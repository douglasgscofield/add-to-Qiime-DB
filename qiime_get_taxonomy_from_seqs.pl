#!/usr/bin/env perl

# Create Qiime-compatible taxonomy from sequence input file
# 
# Modified from:
#
# Script: taxonomy.pl
# Description: Takes an ID file and determines the taxonomy hierarchy for each species 
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
#        sahre001@ucr.edu
# Date: 8.30.11
#       v.1.0  :
#               [x] Make local db
#               [x] generate hierarchy
#               [x] identify non-ranked levels
#               [ ] flag missing levels
#####################################
# Usage: taxonomy.pl [-s] IDfile
#####################################
# -s: flag to show or hide ranks with value "no rank"
#     (include flag to show ranks)
#####################################
# ID file should have the format:
#  someID_someTaxID
# where
#  someID = either Accession or GI
# and
#  someTaxID = NCBI Taxonomy ID
#####################################

use strict;
use warnings;
use Bio::DB::Taxonomy;
use File::Spec;
use Bio::SeqIO;
use Getopt::Long;

my ($name);
my $n_seqs = 0;
my $n_unidentified = 0; # Unidentified organisms
my $n_excluded = 0; # Excluded organisms
my $regexp_excluded; # Regexp to exclude
my $show = 0; # Show or hide ranks with value "no rank"
my $o_qiime = 1; # produce taxon descriptions compatible with qiime
my $o_first = 0; # if there are multiple taxon ids separated by ' ' or ';', take the first rather than failing
my $o_retry = 1; # if there are multiple taxon ids separated by ' ' or ';', try them in turn rather than failing
my $o_taxon = 1;
my $idformat = 'HUDS%0.4u.01FU_%s';  # field 1 is filled with $id1, field 2 is filled with accession
my $id1 = 1;
my $accessionfile;
my $outext = "";
# Handle command-line options
my $debug;

GetOptions( 's|show!' => sub { $show = 1; $outext = '.norank'},
            'q|qiime' => \$o_qiime,
            'noqiime' => sub { $o_qiime = 0 },
            'accessionfile=s' => \$accessionfile,
            'f|first' => \$o_first,
            'idformat' => \$idformat,
            'id1' => \$id1,
            't|taxon' => \$o_taxon,
            'notaxon' => sub { $o_taxon = 0 },
            'retry' => \$o_retry,
            'noretry' => sub { $o_retry = 0 },
            'exclude=s' => \$regexp_excluded,
	        'verbose|debug!' => \$debug,
	    );
die "must use --accessionfile to supply a file of IDs and accession numbers" if $o_qiime and not $accessionfile;
die "must use --idformat to supply an ID format to use" if $o_qiime and not $idformat;
die "only one of --first and --retry may be specified" if $o_first and $o_retry;
my %ACCESSION;
if ($o_qiime) {
    open(my $fh, "<", $accessionfile) or die "could not open $accessionfile: $!";
    while (<$fh>) {
        chomp;
        my @l = split;
        $ACCESSION{$l[0]} = $l[1];
    }
}

my $IDfile = shift @ARGV;


## Function to process a scientific name based on taxonomic rank
#   Returns: the name in angle brackets "<>" if the rank is "no rank"
sub getName {
  my $to = shift;
  my $rank = $to->rank();
  my $res = $to->scientific_name;
  if($rank eq "no rank") {
      $res = join("","<",$to->scientific_name,">");
  }
  return $res;
}

## Make local taxonomy db
my $nodesfile = "./ncbi/nodes.dmp";
my $namefile = "./ncbi/names.dmp";
print STDERR "Loading flatfile taxonomy db from '$nodesfile' and '$namefile'...\n";
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => $nodesfile,
                                   -namesfile => $namefile);
print STDERR "done loading taxonomy db\n";

## Store filename for output files
my @f = split(/\./,$IDfile);
pop(@f);
my $filename = join(".",@f);
$filename = join("",$filename,$outext);

## Get descriptions from ID file
my $seqio = Bio::SeqIO->new(-format => 'fasta',
			    -file   =>  $IDfile);
warn "Creating taxonomy hierarchy..\n" if $debug;

my $outfile_seqs = "$filename\.taxonomy.fa";
my $outfile_taxonomy = "$filename\.taxonomy";
my $outfile_unidentified = "$filename\.unidentified";
my $outfile_excluded = "$filename\.excluded";
my $out = Bio::SeqIO->new(-format=> 'fasta',
			              -file  => ">$outfile_seqs");
open(TAXONOMY,">$outfile_taxonomy");
open(UNIDENTIFIED,">$outfile_unidentified");
print STDERR "sequences to    : $outfile_seqs\n";
print STDERR "taxonomy to     : $outfile_taxonomy\n";
print STDERR "unidentifieds to: $outfile_unidentified\n";
if ($regexp_excluded) {
    open(EXCLUDED,">$outfile_excluded");
    print STDERR "excluded to     : $outfile_excluded    when taxonomic hierarchy matches '$regexp_excluded'\n"
}
# %RANKS keys are acceptable taxonomic ranks found in $curr->rank.  Target of each key is the
# prefix assigned to that rank for the Qiime taxonomy list
sub process_name($) {
    my $n = shift;
    return sprintf($idformat, $id1++, $ACCESSION{$n});
}
my %RANKS = qw/ kingdom k__ phylum p__ class c__ order o__ family f__ genus g__ species s__ /;
while (my $seq = $seqio->next_seq) { 
    printf STDERR "processing sequence $n_seqs\n" if (++$n_seqs % 500) == 0;
    my $display_name = $seq->display_name;
    #print $display_name,"\n";
    ## GI/Accession number is the first item before the first underscore
    ## Organism name/Taxonomy ID is everything after the first underscore
    my ($ID,@descriptors) = split(/\_/, $display_name);
    my $org = join(" ",@descriptors);
    my $taxon_descriptors = $org;
    die "no taxonid found for id '$ID'" if ! @descriptors;
    die "too many taxonid descriptors for id '$ID': $taxon_descriptors" if not ($o_first or $o_retry) and scalar(@descriptors) > 1;
    my @o;
    if ($org =~ /^[\d;]+$/) {
        if ($org =~ /;/) {
            if ($o_first or $o_retry) {
                print STDERR "splitting multi-taxonid descriptor '$org'\n" if $debug;
                @o = split /;/, $org;
                $org = shift @o;
            } else {
                die "taxonid '$org' must be all-numeric; multiple ids may be separated by ';' if --first or --retry specified";
            }
        }
    }

    ## If the input provided is a Taxonomy ID number, 
    ##  get the corresponding organism name
    my $input = $org;
    my $taxonid;
    if ($input =~ /^\d+$/) {  # we need to look up the taxonid to set $input to its scientific_name
        my $tax;
        while ($org) {
            $tax = $taxdb->get_taxon(-taxonid => $org);
            if (defined($tax)) {
                $taxonid = $org;
                $input = $tax->scientific_name;
                last;
            } elsif (! $o_first and @o and $o[0] =~ /^\d+$/) { # more valid taxonids to examine
                print STDERR "taxonid '$org' is undefined, trying '$o[0]'...\n" if $debug;
                $org = shift @o;
            } else { # no valid taxonids
                print STDERR "undefined taxon for taxonid = '$org'\n" if $debug;
                print STDERR "no valid taxonid found in '$taxon_descriptors'\n";
                last;
            }
        } 
    }

    ## Using organism name, generate hierarchy
    my @hierarchy;
    if(my $curr = $taxdb->get_taxon(-name => $input)) {
        my $name = getName($curr);
        #print $name,";";
        #unshift(@hierarchy,$name);

        # could also just check and see if $curr->rank is NULL?
        while($curr) {
            if (exists($RANKS{$curr->rank})) {
                #print "<",$curr->rank,">;";
                $name = getName($curr);
                #print $name,";";
                my $name_no_rank = $name =~ /^</;
                if (not $name_no_rank and $o_qiime) {
                    $name = $RANKS{$curr->rank} . $name;  # add [kpcofgs]__ prefix to rank
                    $name .= " taxonid_".$taxonid if $o_taxon and $taxonid and $curr->rank eq "species";
                }
                ## Default: hide "no rank"
                ## Flag: -s = show "no rank"
                if (not $name_no_rank or $show)  {
                    unshift(@hierarchy,$name) unless $curr->rank eq "species" and not $o_qiime;  # original didn't add species
                }
            }
            $curr = $curr->ancestor;
        }
        my $hierarchy = join(";", @hierarchy);
        if ($regexp_excluded and $hierarchy =~ $regexp_excluded) {
            ++$n_excluded;
            print STDERR "excluding '$ID $input $hierarchy', matches '$regexp_excluded'\n";# if $debug;
            print EXCLUDED "$ID\t$taxon_descriptors\t$hierarchy\n";
            next;
        }
        #print "\n";
        my $final_name = process_name($display_name);
        print TAXONOMY "$final_name\t$hierarchy\n";
        $seq->display_name($final_name);
        $seq->description($hierarchy.($o_qiime ? "" : "."));
        $out->write_seq($seq);      
    } else {
        ## Flag missing taxonomic ranks
        ++$n_unidentified;
        print UNIDENTIFIED "$ID\t$taxon_descriptors\n";
    }
}
print STDERR "Total sequences seen                    : $n_seqs\n";
print STDERR "Total sequences identified and included : ".($id1 - 1)."\n";
print STDERR "Organisms not found in taxonomy database: $n_unidentified".($n_unidentified ? "    See '$outfile_unidentified'." : "")."\n";
if ($regexp_excluded) {
    print STDERR "Organisms excluded by matching '$regexp_excluded': $n_excluded".($n_excluded ? "    See '$outfile_excluded'." : "")."\n";
} else {
    print STDERR "No --exclude option provided.\n";
}

close(UNIDENTIFIED);
close(TAXONOMY);
