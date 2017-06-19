#!/usr/bin/env perl

# Create Qiime-compatible taxonomy from seq id input file, as part of 
#
#     https://github.com/douglasgscofield/add-to-Qiime-DB
#
# Douglas G. Scofield
# Evolutionary Biology Centre, Uppsala University
# douglas.scofield@ebc.uu.se
#
# Heavily modified from Steven Ahrendt's original script found at:
#
#      https://github.com/hyphaltip/mobedac-fungi/blob/master/scripts/taxonomy.pl
#
# which had the original attribution header:
#
## Script: taxonomy.pl
## Description: Takes an ID file and determines the taxonomy hierarchy for each species
## Author: Steven Ahrendt
## email: sahrendt0@gmail.com
##        sahre001@ucr.edu
## Date: 8.30.11
##       v.1.0  :
#
# I've made many modifications to the original script, with the goal of
# adapting the input and output streams for use in producing Qiime-compatible
# output.

use strict;
use warnings;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Taxonomy;
use Bio::SeqIO;

my $n_seqs         = 0;
my $n_taxonomy     = 0;
my $n_output_seqs  = 0;
my $n_unidentified = 0;    # Unidentified organisms
my $n_excluded     = 0;    # Excluded organisms
my @regexp_excluded;       # Regexps to exclude
#push @regexp_excluded, '(^(?!k__Fungi))';
my $n_incomplete         = 0;             # Taxonomy incomplete: includes C or O or F or G and not P
my $n_replaced           = 0;             # Taxonomy replaced from $incompletefile
my $n_redundant          = 0;             # Taxonomically redundantedundant organisms
my $n_truncated          = 0;             # number of sequences truncated for $seq->length > $o_min_to_truncate
my $n_expanded_truncated = 0;             # number of sequences not truncated when $seq->length > $o_min_to_truncate, because of $o_min_after_truncate
my $show                 = 0;             # Show or hide ranks with value "no rank"
my $o_first              = 0;             # if there are multiple taxon ids separated by ' ' or ';', take the first rather than failing
my $o_retry              = 1;             # if there are multiple taxon ids separated by ' ' or ';', try them in turn rather than failing
my $o_taxonid            = 0;             # if taxonid_nnnnn should be appended to the taxonomic hierarchy
my $o_no_redundant       = 0;             # do not produce entries with redundant taxonomy; choose the longest HSP
my $o_db_directory       = "ncbi";
my $o_db_index_directory = "ncbi-indices";
my $taxsep   = '_tax';
my $re_taxsep = qr/$taxsep/;
my $o_usetmp = 0;
my $accessionfile;
my $incompletefile;

# Handle command-line options
my $debug;
my $o_help = 0;

sub usage {
    print STDERR "\n*** Error: @_ \n" if @_;
    print STDERR "
$0: Process blast results and GenBank sequences to produce Qiime-compatible taxonomy and sequence files.
See https://github.com/douglasgscofield/add-to-Qiime-DB for more information.

USAGE:
    $0 [ options ] --accession-file FILE

The FILE given to the --accession-file option may have any name, but must be
produced by the companion script qiime_get_blast_ids_for_genbank.pl.

The NCBI taxonomy DB must be downloaded prior to running this script; see the
Github repository for further information.  For best performance an index should
be generated for it ahead of time, see the Github repository for the procedure.
Locations of the DB and its indices are specified with the --db-directory and
--db-index-directory options.
Command-line options:
  --db-directory DIR         Directory for NCBI taxonomy DB  [$o_db_directory]
  --db-index-directory DIR   Directory for NCBI taxonomy DB indices  [$o_db_index_directory]
  --accession-file FILE      *REQUIRED* File containing blast result columns after 
                             processing with qiime_get_blast_ids_for_genbank.pl
  --incompletefile FILE      NCBI taxonomic hierarchies are often incomplete, missing class,
                             family, etc.  You can replace incomplete hierarchies with complete
                             hierarchies using this file, which should contain the incomplete 
                             hierarchy in column 1, and the complete hierarchy it should be
                             replaced with in column 2.  Note that the hierarchy format should 
                             follow the conventions of this script (hence Qiime).  Each run of
                             the script produces a file suffixed with '.incomplete', containing
                             all incomplete hierarchies.
  --taxonid                  Append the taxonid to the taxonomy string as ' taxonid_xxxx'.  This 
                             is not correct for Qiime input, but can be useful while diagnosing
                             problems and evaluating incomplete taxonomic hierarchies.  [$o_taxonid]
  --redundant                Allow sequences with redundant taxonomic hierarchies.
  --no-redundant             Do not allow sequences with redundant taxonomic hierarchies,
                             pick the longest HSP.  Only one blast hits against sequences with
                             the same taxonomic ID will be kept.  [$o_no_redundant]
  --first                    When multiple taxon descriptors are found, use the 
                             first to find taxonomic information rather than fail.  [$o_first]
  --retry                    When multiple taxon descriptors are found, use each 
                             in turn to find taxonomic information rather than fail.  [$o_retry]
  --no-retry                 When multiple taxon descriptors, fail unless --first specified
  --exclude STRING [ --exclude STRING ... ]
                             Exclude sequences with taxonomic hierarchies that match the 
                             regular expression provided in STRING.  May be specified multiple
                             times.  [ @regexp_excluded ]
  --reset-exclude            Remove all default 'exclude' expressions.  See the default exclude,
                             which removes all but entries for Kingdom Fungi, for an example of
                             a Perl regexp that excludes all but a specific kingdom.
  --usetmp                   Include a temporary identifier (process ID) in filenames [$o_usetmp]
  --no-usetmp                Do not include a temporary identifier ...
  --verbose | --debug        Produce more message output than you probably care to see.
";
    exit 1;
}

GetOptions(
    'db-directory=s'       => \$o_db_directory,
    'db-index-directory=s' => \$o_db_index_directory,
    'accession-file=s'     => \$accessionfile,
    'incompletefile=s'     => \$incompletefile,
    'taxonid'              => \$o_taxonid,
    'no-taxonid'           => sub { $o_taxonid = 0 },
    'redundant'            => sub { $o_no_redundant = 0 },
    'no-redundant'         => \$o_no_redundant,
    'first'                => \$o_first,
    'retry'                => \$o_retry,
    'no-retry'             => sub { $o_retry = 0 },
    'reset-exclude'        => sub { @regexp_excluded = () },
    'exclude=s'            => \@regexp_excluded,
    'usetmp'               => \$o_usetmp,
    'no-usetmp'            => sub { $o_usetmp = 0 },
    'verbose|debug!'       => \$debug,
    'help|?'               => \$o_help,
);
usage() if $o_help;
usage("must provide --accession-file only")
  if not $accessionfile or @ARGV != 0;
usage("only one of --first and --retry may be specified")     if $o_first and $o_retry;

my %ACCESSION;

sub force_orientation($$) {
    my ( $s, $e ) = @_;
    return $s > $e ? ( $e, $s ) : ( $s, $e );    # swap if on reverse strand
}

# Read --incompletefile file containing replacement info for incomplete hierarchies

my %INCOMPLETE;
my $n_incompletefile_lines = 0;
if ($incompletefile) {
    open( my $fh, "<", $incompletefile ) or die "could not open $incompletefile $!";
    while (<$fh>) {
        ++$n_incompletefile_lines;
        chomp;
        my @l = split /\t/;
        $INCOMPLETE{ $l[0] } = $l[1];
    }
    print STDERR "\nLoaded incomplete taxonomic hierarchies from $incompletefile containing $n_incompletefile_lines lines\n\n";
}

# Read local taxonomic database

my $nodesfile = "$o_db_directory/nodes.dmp";
my $namefile  = "$o_db_directory/names.dmp";
print STDERR "\nLoading flatfile taxonomy db from '$nodesfile' and '$namefile', with indices in '$o_db_index_directory/' ...\n";
my $taxdb = Bio::DB::Taxonomy->new(
    -source    => 'flatfile',
    -directory => $o_db_index_directory,    # the location of the index files
    -nodesfile => $nodesfile,
    -namesfile => $namefile
);
print STDERR "done loading taxonomy db\n\n";

print STDERR "Calculating taxonomy using just the IDs in $accessionfile\n\n";

#
# Utility functions for taxonomic info

## Function to process a scientific name based on taxonomic rank
#   Returns: the name in angle brackets "<>" if the rank is "no rank"
sub getName {
    my $to   = shift;
    my $rank = $to->rank();
    my $res  = $to->scientific_name;
    if ( $rank eq "no rank" ) {
        $res = join( "", "<", $to->scientific_name, ">" );
    }
    return $res;
}

# %RANKS keys are acceptable taxonomic ranks found in $curr->rank.  Target of each key is the
# prefix assigned to that rank for the Qiime taxonomy list

my %RANKS = qw/ kingdom k__ phylum p__ class c__ order o__ family f__ genus g__ species s__ /;

sub find_taxonid($) {
    my $seq = shift;
    my ( $ID, @descriptors ) = split $re_taxsep, $seq->display_name;
    return process_taxon_descriptors($ID, @descriptors);
}

sub process_taxon_descriptors($@) {
    my ($ID, @descriptors) = @_;
    my $org = join( " ", @descriptors );
    my $taxon_descriptors = $org;
    die "no taxonid found for id '$ID'" if !@descriptors;
    die "too many taxonid descriptors for id '$ID': $taxon_descriptors"
      if not( $o_first or $o_retry )
      and scalar(@descriptors) > 1;
    my @o;
    if ( $org =~ /^[\d;]+$/ ) {

        if ( $org =~ /;/ ) {
            if ( $o_first or $o_retry ) {
                print STDERR "splitting multi-taxonid descriptor '$org'\n" if $debug;
                @o = split /;/, $org;
                $org = shift @o;
            }
            else {
                die "taxonid '$org' must be all-numeric; multiple ids may be separated by ';' if --first or --retry specified";
            }
        }
    }

    # If the input provided is a Taxonomy ID number, get the corresponding organism name
    my $input   = $org;
    my $taxonid = "";
    if ( $input =~ /^\d+$/ ) {    # we need to look up the taxonid to set $input to its scientific_name
        my $tax;
        while ($org) {
            $tax = $taxdb->get_taxon( -taxonid => $org );
            if ( defined($tax) ) {
                $taxonid = $org;
                $input   = $tax->scientific_name;
                last;
            }
            elsif ( !$o_first and @o and $o[0] =~ /^\d+$/ ) {    # more valid taxonids to examine
                print STDERR "taxonid '$org' is undefined, trying '$o[0]'...\n" if $debug;
                $org = shift @o;
            }
            else {                                               # no valid taxonids
                print STDERR "undefined taxon for taxonid = '$org'\n" if $debug;
                print STDERR "no valid taxonid found in '$taxon_descriptors'\n";
                last;
            }
        }
    }
    return ( $ID, $input, $taxonid, $taxon_descriptors );
}

sub generate_hierarchy($$$) {
    my $curr       = shift;
    my $taxonid    = shift;
    my $ranks_seen = shift;
    return ""                                                       if not $curr;         # undefined taxonomy
    die "generate_hierarchy() requires valid taxonid for arg 2"     if not $taxonid;
    die "generate_hierarchy() requires reference to hash for arg 3" if not $ranks_seen;
    my @hierarchy;
    %{$ranks_seen} = map { $_ => 0 } keys %RANKS;

    while ($curr) {
        if ( exists $RANKS{ $curr->rank } ) {
            my $name         = getName($curr);
            my $name_no_rank = $name =~ /^</;
            if ( not $name_no_rank ) {
                $name = $RANKS{ $curr->rank } . $name;                                    # add [kpcofgs]__ prefix to rank
                $ranks_seen->{ $curr->rank }++;
                $name .= " taxonid_" . $taxonid if $taxonid and $curr->rank eq "species";
                unshift @hierarchy, $name;
            }
        }
        $curr = $curr->ancestor;
    }
    $ranks_seen->{is_incomplete}  = 0;
    $ranks_seen->{was_incomplete} = 0;
    if (   !$ranks_seen->{kingdom}
        or !$ranks_seen->{phylum}
        or !$ranks_seen->{class}
        or !$ranks_seen->{order}
        or !$ranks_seen->{family}
        or !$ranks_seen->{genus} )
    {
        $ranks_seen->{is_incomplete} = 1;
    }
    my $hierarchy = join( ";", @hierarchy );
    if ( exists $INCOMPLETE{$hierarchy} ) {
        $hierarchy                    = $INCOMPLETE{$hierarchy};
        $ranks_seen->{was_incomplete} = 1;
        $ranks_seen->{is_incomplete}  = 0;
    }

    return $hierarchy;
}

## Store filename for output files
#my @f = split( /\./, $IDfile );
my @f = split( /\./, $accessionfile );
pop(@f);
my $filename = join( ".", @f );
$filename = join( "", $filename, ( $o_usetmp ? ( "." . $$ ) : "" ) );

# Get descriptions from ID file

my $outfile_taxonomy      = "$filename.taxonomy";
my $outfile_unidentified  = "$filename.unidentified";
my $outfile_redundant     = "$filename.redundant";
my $outfile_excluded      = "$filename.excluded";
my $outfile_incomplete    = "$filename.incomplete";
my $outfile_wasincomplete = "$filename.wasincomplete";

open( TAXONOMY,     ">$outfile_taxonomy" );
open( UNIDENTIFIED, ">$outfile_unidentified" );
open( INCOMPLETE,   ">$outfile_incomplete" );
open( REDUNDANT,    ">$outfile_redundant" );
open( EXCLUDED,     ">$outfile_excluded" );

print STDERR "
accessions in file: $accessionfile
taxonomy          : $outfile_taxonomy
unidentified      : $outfile_unidentified
incomplete        : $outfile_incomplete
redundant         : $outfile_redundant
excluded          : $outfile_excluded when matching /" . join( "/ /", @regexp_excluded ) . "/
";

# Read accession file containing info in HSPs of blast hits

# Go through accessions, find taxonid, calculate hierarchy, add hierarchy info to
# ACCESSION, and add to ENTRIES list for each hierarchy observed

print STDERR "\nPass 1: Read accessions, calculate taxonomic hierarchies, note redundancies and exclusions\n\n";

my %ENTRIES;

my $n_accessionfile_lines = 0;
{
    open( my $fh, "<", $accessionfile ) or die "could not open $accessionfile: $!";
    while (<$fh>) {
        ++$n_accessionfile_lines;
        chomp;
        my @l          = split /\t/;
        my $accv_taxonid = $l[0] . $taxsep . $l[1];
        my ( $start, $end ) = force_orientation( $l[2], $l[3] );
        $ACCESSION{ $accv_taxonid }{key}          = $accv_taxonid;
        $ACCESSION{ $accv_taxonid }{accession}    = $l[0];
        $ACCESSION{ $accv_taxonid }{taxonid}      = $l[1];
        $ACCESSION{ $accv_taxonid }{start}        = $start;
        $ACCESSION{ $accv_taxonid }{end}          = $end;
        $ACCESSION{ $accv_taxonid }{hsplen}       = $end - $start + 1;
        $ACCESSION{ $accv_taxonid }{is_redundant} = 0;

        my ( $ID, $input, $taxonid, $taxon_descriptors ) = process_taxon_descriptors($l[0], $l[1]);
        my $curr = $taxdb->get_taxon( -name => $input );
        my %ranks_seen;
        my $hierarchy = generate_hierarchy( $curr, $taxonid, \%ranks_seen );
        print STDERR "ENTRIES: for accv_taxonid " . $accv_taxonid . " hierarchy '$hierarchy'\n"
          if $debug;

        $ACCESSION{ $accv_taxonid }{hierarchy}         = $hierarchy;
        $ACCESSION{ $accv_taxonid }{is_incomplete}     = $ranks_seen{is_incomplete};
        $ACCESSION{ $accv_taxonid }{was_incomplete}    = $ranks_seen{was_incomplete};
        $ACCESSION{ $accv_taxonid }{input}             = $input;
        $ACCESSION{ $accv_taxonid }{ID}                = $ID;
        $ACCESSION{ $accv_taxonid }{taxonid}           = $taxonid;
        $ACCESSION{ $accv_taxonid }{taxon_descriptors} = $taxon_descriptors;

        if ( exists $ENTRIES{$hierarchy} ) {
            print STDERR "ENTRIES: hierarchy '$hierarchy' redundancy found for accv_taxonid " . $accv_taxonid . "\n"
              if $debug;
            push @{ $ENTRIES{$hierarchy} }, $ACCESSION{ $accv_taxonid };
        }
        else {
            print STDERR "ENTRIES: hierarchy '$hierarchy' created for accv_taxonid " . $accv_taxonid . "\n"
              if $debug;
            $ENTRIES{$hierarchy} = [ $ACCESSION{ $accv_taxonid } ];
        }

        # Exclude by matching sequences given in --exclude
        if (@regexp_excluded) {
            $ACCESSION{ $accv_taxonid }{excluded} = 0;
            foreach my $regexp (@regexp_excluded) {
                if ( $hierarchy =~ $regexp ) {
                    $ACCESSION{ $accv_taxonid }{excluded} = $regexp;
                    last;
                }
            }
        }
    }
}
print STDERR "\nLoaded accessions from $accessionfile containing $n_accessionfile_lines lines\n\n";

print STDERR "\nChecking for redundant taxa\n\n\%ENTRIES contains " . scalar( keys %ENTRIES ) . " keys\n\n";
foreach ( keys %ENTRIES ) {
    my $n = scalar( @{ $ENTRIES{$_} } );
    next if $n == 1;                        # only one entry with this hierarchy, leave it alone
    die "ENTRIES{$_} length is 0" if !$n;

    print STDERR "\$ENTRIES{$_} has redundancy: $n entries" . ($o_no_redundant ? " so keeping largest HSP" : " but keeping all entries") . "\n";

    my ( $maxlen, $maxref ) = ( 0, 0 );
    my $e = $_;
    foreach my $a ( @{ $ENTRIES{$e} } ) {
        if ( $a->{hsplen} > $maxlen ) {
            print STDERR "    " . $a->{key} . " longest hsp " . $a->{hsplen} . ", old " . ( $maxref ? $maxref->{key} : "<init>" ) . " $maxlen\n";
            $maxref->{is_redundant} = 1 if $maxref;    # invalidate the old larger
            $maxlen                 = $a->{hsplen};
            $maxref                 = $a;
            $maxref->{is_redundant} = 0;               # validate the new larger
        }
        else {
            $a->{is_redundant} = 1;                    # invalidate the smaller
        }
    }
}

# reopen accession file and generate new taxonomy after excluding redundants etc.

print STDERR "\nPass 2: Read accessions again, output those not redundant or excluded, truncate where necessary\n\n";

{
    open( my $fh, "<", $accessionfile ) or die "could not open $accessionfile: $!";
    while (<$fh>) {
        chomp;
        my @l          = split /\t/;
        my $accv_taxonid = $l[0] . $taxsep . $l[1];

        die "cannot find accession for ", $accv_taxonid
          if not exists $ACCESSION{ $accv_taxonid };

        my $acc = $ACCESSION{ $accv_taxonid };

        my ( $hierarchy, $ID, $input, $taxonid, $taxon_descriptors ) =
          @{$acc}{qw/hierarchy ID input taxonid taxon_descriptors/};

        # If we could not detect a taxonomic hierarchy earlier...
        if ( not $acc->{hierarchy} ) {    # this is defined to be ""
            ++$n_unidentified;
            print UNIDENTIFIED "$ID\t$taxon_descriptors\n";
            next;
        }

        # Has this accession already been declared redundant?
        die "\%ENTRIES value does not exist for hierarchy '$hierarchy'"
          if not exists $ENTRIES{$hierarchy};
        if ( $acc->{is_redundant} ) {
            ++$n_redundant;
            print REDUNDANT "ACCESSION_REDUNDANT\t$accv_taxonid\t$hierarchy\n";
            next if $o_no_redundant;
        }

        # Exclude by matching pattern given in --exclude
        if ( $acc->{excluded} ) {
            ++$n_excluded;
            print STDERR "HIERARCHY_EXCLUDED '$ID', '$input', '$hierarchy': matches '" . $acc->{excluded} . "'\n";    # if $debug;
            print EXCLUDED "HIERARCHY_EXCLUDED\t$ID\t$taxon_descriptors\t$hierarchy\n";
            next;
        }

        # We expect to output an accession with this hierarchy

        # This is the name we will give the sequence for the final output
        my $final_name = $accv_taxonid;

        # If we modified the taxonomy, or it was incomplete, report that
        if ( $acc->{was_incomplete} ) {
            ++$n_replaced;    # don't set in generate_hierarchy(), that might be used at other times
            print INCOMPLETE "HIERARCHY_REPLACED\t$final_name\t$hierarchy\n";
        }
        elsif ( $acc->{is_incomplete} ) {

            # incomplete taxonomy, note it
            ++$n_incomplete;
            print INCOMPLETE "HIERARCHY_INCOMPLETE\t$final_name\t$hierarchy\n";
        }

        # remove the taxonid, it helped us with redundancy but the user doesn't want it on output
        if ( not $o_taxonid ) {
            $hierarchy =~ s/ taxonid_\d+$//;
        }

        ++$n_taxonomy;
        print TAXONOMY "$final_name\t$hierarchy\n";

    } # read accession line
} # code block

printf STDERR "
%8d accessions in file: $accessionfile
%8d taxonomy          : $outfile_taxonomy
%8d unidentified      : $outfile_unidentified
%8d incomplete        : $outfile_incomplete with 'HIERARCHY_INCOMPLETE'
%8d replaced          : $outfile_incomplete with 'HIERARCHY_REPLACED'
%8d redundant         : $outfile_redundant
%8d excluded          : $outfile_excluded when matching /" . join( "/ /", @regexp_excluded ) . "/
", $n_accessionfile_lines, $n_taxonomy, $n_unidentified, $n_incomplete, $n_replaced, $n_redundant, $n_excluded;

