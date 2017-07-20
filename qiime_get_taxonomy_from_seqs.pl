#!/usr/bin/env perl

# Create Qiime-compatible taxonomy from sequence input file, as part of 
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
push @regexp_excluded, '(^(?!k__Fungi))';
my $n_incomplete         = 0;             # Taxonomy incomplete: includes C or O or F or G and not P
my $n_replaced           = 0;             # Taxonomy replaced from $incompletefile
my $n_redundant          = 0;             # Taxonomically redundantedundant organisms
my $n_truncated          = 0;             # number of sequences truncated for $seq->length > $o_min_to_truncate
my $n_expanded_truncated = 0;             # number of sequences not truncated when $seq->length > $o_min_to_truncate, because of $o_min_after_truncate
my $show                 = 0;             # Show or hide ranks with value "no rank"
my $o_first              = 0;             # if there are multiple taxon ids separated by ' ' or ';', take the first rather than failing
my $o_retry              = 1;             # if there are multiple taxon ids separated by ' ' or ';', try them in turn rather than failing
my $o_taxonid            = 0;             # if taxonid_nnnnn should be appended to the taxonomic hierarchy
my $o_no_redundant       = 1;             # do not produce entries with redundant taxonomy; choose the longest HSP
my $o_min_to_truncate    = 400;           # if the GenBank sequence is longer than this, truncate it to the extent of the blast HSP
my $o_min_after_truncate = 300;           # if the GenBank sequence would be shorter than this after truncation, expand it until it is this long
my $o_db_directory       = "ncbi";
my $o_db_index_directory = "ncbi-indices";
my $taxsep   = '_tax';
my $re_taxsep = qr/$taxsep/;
my $idformat = 'ID_%0.4u.01FU_%s';       # field 1 is filled with $id1, field 2 is filled with accession
my $id1      = 1;
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
    $0 [ options ] --accession-file FILE FILE.fa

The FILE given to the --accession-file option may have any name, but must be
produced by the companion script qiime_get_blast_ids_for_genbank.pl.
'FILE.fa' may have any name, but must be a Fasta-format file produced by the
companion script qiime_get_genbank_seqs.pl.

The NCBI taxonomy DB must be downloaded prior to running this script, and an
index should be generated; see the Github repository for further information.
Locations of the DB and its indices are specified with the --db-directory and
--db-index-directory options, if necessary.

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
  --min-to-truncate INT      Sometimes sequences representing the best hit for a query within 
                             GenBank are very large; using this option specifies the minimum
                             length (bp) a GenBank sequence must be to be truncated in length to the
                             extent of the HSP of the original blast hit.  The region surrounding
                             the HSP can be expanded after truncation using --min-after-truncate.
                             [$o_min_to_truncate]
  --min-after-truncate INT   Minimum length (bp) of region around blast HSP, grown equally 
                             around the HSP if necessary.  [$o_min_after_truncate]
  --taxonid                  Append the taxonid to the taxonomy string as ' taxonid_xxxx'.  This 
                             is not correct for Qiime input, but can be useful while diagnosing
                             problems and evaluating incomplete taxonomic hierarchies.  [$o_taxonid]
  --redundant                Allow sequences with redundant taxonomic hierarchies.
  --no-redundant             Do not allow sequences with redundant taxonomic hierarchies, 
                             pick the longest HSP.  Only one blast hits against sequences with 
                             the same taxonomic ID will be kept.  [$o_no_redundant]
  --first                    When multiple taxon descriptors are found for a sequence, use the 
                             first to find taxonomic information rather than fail.  [$o_first]
  --retry                    When multiple taxon descriptors are found for a sequence, use each 
                             in turn to find taxonomic information rather than fail.  [$o_retry]
  --no-retry                 When multiple taxon descriptors, fail unless --first specified
  --idformat STRING          Output sequence IDs with this format.  The format should contain
                             two fields, an unsigned integer filled with sequential numbers
                             starting with --id1, e.g. '%0.4u' and a string field '%s' filled
                             with the GenBank accession.  [$idformat]
  --id1 INT                  Begin numbering sequences with this number, see --idformat.  [$id1]
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
    'min-to-truncate=i'    => \$o_min_to_truncate,
    'min-after-truncate=i' => \$o_min_after_truncate,
    'taxonid'              => \$o_taxonid,
    'no-taxonid'           => sub { $o_taxonid = 0 },
    'redundant'            => sub { $o_no_redundant = 0 },
    'no-redundant'         => \$o_no_redundant,
    'first'                => \$o_first,
    'retry'                => \$o_retry,
    'no-retry'             => sub { $o_retry = 0 },
    'idformat=s'           => \$idformat,
    'id1=i'                => \$id1,
    'reset-exclude'        => sub { @regexp_excluded = () },
    'exclude=s'            => \@regexp_excluded,
    'usetmp'               => \$o_usetmp,
    'no-usetmp'            => sub { $o_usetmp = 0 },
    'verbose|debug!'       => \$debug,
    'help|?'               => \$o_help,
);

usage() if $o_help;
usage("must provide --accession-file and one file of sequences") if not $accessionfile or @ARGV != 1;
usage("must use --idformat to supply an ID format to use")    if not $idformat;
usage("only one of --first and --retry may be specified")     if $o_first and $o_retry;
usage("--id1 must be greater than or equal to 0")             if $id1 < 0;
usage("--min-to-truncate must be greater than or equal to 0") if $o_min_to_truncate < 0;
usage("--min-after-truncate must be less than or equal to --min-to-truncate") if $o_min_after_truncate > $o_min_to_truncate;

# Read accession file containing info in HSPs of blast hits

my %ACCESSION;

sub force_orientation($$) {
    my ( $s, $e ) = @_;
    return $s > $e ? ( $e, $s ) : ( $s, $e );    # swap if on reverse strand
}

open( my $fh, "<", $accessionfile ) or die "could not open $accessionfile: $!";
my $n_accessionfile_lines = 0;
while (<$fh>) {
    chomp;
    ++$n_accessionfile_lines;
    my @l          = split /\t/;
    my $accv_taxonid = $l[0] . $taxsep . $l[1];
    my ( $start, $end ) = force_orientation( $l[2], $l[3] );
    $ACCESSION{$accv_taxonid}{key}          = $accv_taxonid;
    $ACCESSION{$accv_taxonid}{accession}    = $l[0];
    $ACCESSION{$accv_taxonid}{taxonid}      = $l[1];
    $ACCESSION{$accv_taxonid}{start}        = $start;
    $ACCESSION{$accv_taxonid}{end}          = $end;
    $ACCESSION{$accv_taxonid}{hsplen}       = $end - $start + 1;
    $ACCESSION{$accv_taxonid}{is_redundant} = 0;
}
print STDERR "\nLoaded accessions from $accessionfile containing $n_accessionfile_lines lines\n\n";

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

sub process_name($) {
    my $n = shift;
    return sprintf( $idformat, $id1++, $ACCESSION{$n}{accession} );
}

# %RANKS keys are acceptable taxonomic ranks found in $curr->rank.  Target of each key is the
# prefix assigned to that rank for the Qiime taxonomy list

my %RANKS = qw/ kingdom k__ phylum p__ class c__ order o__ family f__ genus g__ species s__ /;

sub find_taxonid($) {
    my $seq = shift;
    my ( $ID, @descriptors ) = split $re_taxsep, $seq->display_name;
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

my $IDfile = shift @ARGV;
-f $IDfile or die "file of sequences with IDs '$IDfile' does not exist";

## Store filename for output files
my @f = split( /\./, $IDfile );
pop(@f);
my $filename = join( ".", @f );
$filename = join( "", $filename, ( $o_usetmp ? ( "." . $$ ) : "" ) );

# Get descriptions from ID file

my $outfile_seqs          = "$filename.taxonomy.fa";
my $outfile_taxonomy      = "$filename.taxonomy";
my $outfile_unidentified  = "$filename.unidentified";
my $outfile_redundant     = "$filename.redundant";
my $outfile_excluded      = "$filename.excluded";
my $outfile_incomplete    = "$filename.incomplete";
my $outfile_wasincomplete = "$filename.wasincomplete";
my $outfile_truncated     = "$filename.truncated";
my $out                   = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => ">$outfile_seqs"
);
open( TAXONOMY,     ">$outfile_taxonomy" );
open( UNIDENTIFIED, ">$outfile_unidentified" );
open( INCOMPLETE,   ">$outfile_incomplete" );
open( REDUNDANT,    ">$outfile_redundant" );
open( EXCLUDED,     ">$outfile_excluded" );
open( TRUNCATED,    ">$outfile_truncated" );

print STDERR "
input sequences   : $IDfile
output sequences  : $outfile_seqs
taxonomy          : $outfile_taxonomy
unidentified      : $outfile_unidentified
incomplete        : $outfile_incomplete
redundant         : $outfile_redundant
excluded          : $outfile_excluded when matching /" . join( "/ /", @regexp_excluded ) . "/
>$o_min_to_truncate bp truncated : $outfile_truncated
...then expanded to $o_min_after_truncate bp if possible
";

# Go through fasta sequences, find taxonid, calculate hierarchy, add hierarchy info to
# ACCESSION, and add to ENTRIES list for each hierarchy observed

print STDERR "\nPass 1: Read sequences, calculate taxonomic hierarchies, note redundancies and exclusions\n\n";

my %ENTRIES;

my $seqio = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => $IDfile
);
while ( my $seq = $seqio->next_seq ) {
    printf STDERR "Pass 1: processing sequence $n_seqs ...\n" if ( ++$n_seqs % 500 ) == 0;
    my ( $ID, $input, $taxonid, $taxon_descriptors ) = find_taxonid($seq);
    my $curr = $taxdb->get_taxon( -name => $input );
    my %ranks_seen;
    my $hierarchy = generate_hierarchy( $curr, $taxonid, \%ranks_seen );
    print STDERR "ENTRIES: for sequence " . $seq->display_name . " hierarchy '$hierarchy'\n"
      if $debug;

    # verify %ACCESSION entry
    if ( $ACCESSION{ $seq->display_name }{key} ne $seq->display_name ) {
        print STDERR "mismatch between expected key '$seq->display_name' and \$ACCESSION{\$seq->display_name}{key}, Dumper() shows:\n";
        print STDERR Dumper( \$ACCESSION{ $seq->display_name } );
    }
    $ACCESSION{ $seq->display_name }{hierarchy}         = $hierarchy;
    $ACCESSION{ $seq->display_name }{is_incomplete}     = $ranks_seen{is_incomplete};
    $ACCESSION{ $seq->display_name }{was_incomplete}    = $ranks_seen{was_incomplete};
    $ACCESSION{ $seq->display_name }{input}             = $input;
    $ACCESSION{ $seq->display_name }{ID}                = $ID;
    $ACCESSION{ $seq->display_name }{taxonid}           = $taxonid;
    $ACCESSION{ $seq->display_name }{taxon_descriptors} = $taxon_descriptors;
    if ( exists $ENTRIES{$hierarchy} ) {
        print STDERR "ENTRIES: hierarchy '$hierarchy' redundancy found for sequence " . $seq->display_name . "\n"
          if $debug;
        push @{ $ENTRIES{$hierarchy} }, $ACCESSION{ $seq->display_name };
    }
    else {
        print STDERR "ENTRIES: hierarchy '$hierarchy' created for sequence " . $seq->display_name . "\n"
          if $debug;
        $ENTRIES{$hierarchy} = [ $ACCESSION{ $seq->display_name } ];
    }

    # Exclude by matching sequences given in --exclude
    if (@regexp_excluded) {
        $ACCESSION{ $seq->display_name }{excluded} = 0;
        foreach my $regexp (@regexp_excluded) {
            if ( $hierarchy =~ $regexp ) {
                $ACCESSION{ $seq->display_name }{excluded} = $regexp;
                last;
            }
        }
    }
}
$seqio->close();

# go through ENTRIES, and all but one of the sequence entries (the longest HSP) as redundant

print STDERR "\nChecking for redundant taxa\n\n\%ENTRIES contains " . scalar( keys %ENTRIES ) . " keys\n\n";
foreach ( keys %ENTRIES ) {
    my $n = scalar( @{ $ENTRIES{$_} } );
    next if $n == 1;                        # only one entry with this hierarchy, leave it alone
    die "ENTRIES{$_} length is 0" if !$n;

    print STDERR "\$ENTRIES{$_} has redundancy: $n entries\n";

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

# reopen sequence file and generate new sequences and taxonomy after excluding redundants etc.

print STDERR "\nPass 2: Read sequences again, output those not redundant or excluded, truncate where necessary\n\n";

$seqio = Bio::SeqIO->new(
    -format => 'fasta',
    -file   => $IDfile
);

$n_seqs = 0;
while ( my $seq = $seqio->next_seq ) {

    printf STDERR "Pass 2: processing sequence $n_seqs ...\n" if ( ++$n_seqs % 500 ) == 0;

    die "cannot find accession for ", $seq->display_name
      if not exists $ACCESSION{ $seq->display_name };

    my $acc = $ACCESSION{ $seq->display_name };

    my ( $hierarchy, $ID, $input, $taxonid, $taxon_descriptors ) =
      @{$acc}{qw/hierarchy ID input taxonid taxon_descriptors/};

    # If we could not detect a taxonomic hierarchy earlier...
    if ( not $acc->{hierarchy} ) {    # this is defined to be ""
        ++$n_unidentified;
        print UNIDENTIFIED "$ID\t$taxon_descriptors\n";
        next;
    }

    # Has this sequence already been declared redundant?
    die "\%ENTRIES value does not exist for hierarchy '$hierarchy'"
      if not exists $ENTRIES{$hierarchy};
    if ( $acc->{is_redundant} ) {
        ++$n_redundant;
        print STDERR "SEQUENCE_REDUNDANT:\t" . $seq->display_name . "\t" . $seq->length . "\t$hierarchy\n"
          if $debug;
        print REDUNDANT "SEQUENCE_REDUNDANT\t" . $seq->display_name . "\t" . $seq->length . "\t$hierarchy\n";
        next;
    }

    # Exclude by matching sequences given in --exclude
    if ( $acc->{excluded} ) {
        ++$n_excluded;
        print STDERR "HIERARCHY_EXCLUDED '$ID', '$input', '$hierarchy': matches '" . $acc->{excluded} . "'\n";    # if $debug;
        print EXCLUDED "HIERARCHY_EXCLUDED\t$ID\t$taxon_descriptors\t$hierarchy\n";
        next;
    }

    # We expect to output a sequence with this hierarchy

    # This is the name we will give the sequence for the final output
    my $final_name = process_name( $seq->display_name );

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

    # Truncate subject sequence if necessary
    if ( $o_min_to_truncate and $seq->length > $o_min_to_truncate ) {
        ++$n_truncated;
        my $slen = $seq->length;
        my ( $hspstart, $hspend ) = @{ $ACCESSION{ $seq->display_name } }{ 'start', 'end' };
        my $hsplen = $hspend - $hspstart + 1;
        my ( $start, $end ) = ( $hspstart, $hspend );
        my $d      = $o_min_after_truncate - $hsplen;
        my $status = "as_hsp";
        if ( $d > 0 ) {    # we are short of --min-after-truncate, expand the truncation site evenly both sides
            ++$n_expanded_truncated;
            $status = "expanded";
            ++$d if $d % 2;    # make even for cleaner math
            $start = ( $start <= $d / 2 ) ? 1 : ( $start - ( $d / 2 ) );
            my $desired_end = $start + $o_min_after_truncate - 1;
            $end = ( $desired_end > $slen ) ? $slen : $desired_end;
        }
        $seq->seq( $seq->subseq( $start, $end ) );
        my $newlen = $seq->length;
        print STDERR "Truncating " . $seq->display_name . " from $slen to $newlen ($start-$end), HSP is $hsplen ($hspstart-$hspend)\t$status\n"
          if $debug;
        print TRUNCATED "TRUNCATED\t$final_name\t"
          . $seq->display_name
          . "\torig:$slen\thsp:$hsplen:$hspstart-$hspend\t$status\tnew:$newlen\tfrom:$start-$end\t$status\n";
    }

    # remove the taxonid, it helped us with redundancy but the user doesn't want it on output
    if ( not $o_taxonid ) {
        $hierarchy =~ s/ taxonid_\d+$//;
    }

    ++$n_taxonomy;
    print TAXONOMY "$final_name\t$hierarchy\n";

    $seq->display_name($final_name);
    $seq->description($hierarchy);
    ++$n_output_seqs;
    $out->write_seq($seq);

}

printf STDERR "
%8d input sequences   : $IDfile
%8d output sequences  : $outfile_seqs
%8d taxonomy          : $outfile_taxonomy
%8d unidentified      : $outfile_unidentified
%8d incomplete        : $outfile_incomplete with 'HIERARCHY_INCOMPLETE'
%8d replaced          : $outfile_incomplete with 'HIERARCHY_REPLACED'
%8d redundant         : $outfile_redundant
%8d excluded          : $outfile_excluded when matching /" . join( "/ /", @regexp_excluded ) . "/
%8d >$o_min_to_truncate bp truncated : $outfile_truncated
%8d ...then expanded to $o_min_after_truncate bp if possible
", $n_seqs, $n_output_seqs, $n_taxonomy, $n_unidentified, $n_incomplete, $n_replaced, $n_redundant, $n_excluded, $n_truncated, $n_expanded_truncated;

