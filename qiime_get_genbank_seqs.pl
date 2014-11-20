#!/usr/bin/env perl 

    eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell


=head1 NAME

qiime_get_genbank_seqs - Retrieve Genbank sequences and produce output useful
                         for augmenting a Qiime reference database.
 
This script is derived from the BioPerl script bp_download_query_genbank.pl.
Differences include restricting to FASTA format and restricting to the
nucleotide database, and the additional options --accession/--noaccession,
--taxonomy/--notaxonomy, and --primary for controlling output and naming of
output sequences, removal of the query file option, and removal of the output 
file option.

Two files are output: a file of FASTA-format sequences to STDOUT, and a file of
IDs to STDERR.

The IDs provided for the query can include additional taxonomy information tacked
onto the end (see the --taxonomy option) and can include additional information
after the IDs separated by whitespace.  The ID for query is taken from just the
first entry (possibly modified according to --taxonomy).  The IDs output to stderr
are the full ID provided, including any information following whitespace, in the
same order as the sequences in the sequence file written to STDOUT.  Note that
because of the nature of Genbank queries, the output is not likely to be in the
same order as the query IDs in the gifile, but the files written to STDOUT and
STDERR *will be* in the same order.

=head1 USAGE

qiime_get_genbank_seqs --query 18481572 > seqs.fa 2> seqs_ids.txt

qiime_get_genbank_seqs --gifile file_with_gis.txt > seqs.fa 2>seqs_ids.txt

=head2 Other options

Provide ONE of:

=over

=item -q | --query string

string is an ID to use to query GenBank, or

=back

=over

=item --gi --gis --gifile file 

where file contains a list of GIs to download.

=back

This script is modified to support its use with augmenting a Qiime reference database.
The --taxonomy and --accession options are On by default for this reason.

=over

=item -t --taxonomy

query is in the form of gi_taxonid, as also expected by taxonomy.pl the sequence display ID will also be gi_taxonid.  Only the first whitespace-separated field will be used

=item --notaxonomy

used to turn off --taxonomy option, which is On by default

=item -a --accession

write table of ID <TAB> accession number to STDERR

=item --noaccession

used to turn off --accession option, which is On by default

=item -p --primary

prefix the primary ID to the display_id before writing the sequence CANNOT co-occur with --taxonomy if it co-occurs with --verbose, STDERR will be all mixed up

=item -v --verbose

debugging output

=back
 
=head2 Query options

=over

=item --maxids maximum number of IDs to retrieve in a set (default 100)

=item --reldate 

=item --maxdate maxdate for a record

=item --mindate minimum date for record

=item --datetype edat or mdat (entered or modified)

=back

=head1 AUTHOR Jason Stajich

Jason Stajich, jason-AT-bioperl.org

Modified by Douglas Scofield to do Qiime-oriented work, douglas.scofield@ebc.uu.se

=cut

use strict;
use warnings;
use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::DB::Query::GenBank;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my ($format,$debug,$primary,$accession,$taxonomy,%options);
# for Qiime work, all four of these should remain unchanged
$accession = 1;
$taxonomy = 1;
$format = 'fasta';
$options{'-db'} = 'nucleotide';

$options{'-maxids'} = '100';
my $gifile;
GetOptions(
		      'h|help' => sub { exec('perldoc', $0); 
									exit(0);
								},
			  'v|verbose'       => \$debug,
			  'gi|gifile|gis:s' => \$gifile,
			  'p|primary'       => \$primary,
			  'a|accession'     => \$accession,
              'noaccession'     => sub { $accession = 0 },
              't|taxonomy'      => \$taxonomy,
              'notaxonomy'      => sub { $taxonomy = 0 },
			  # DB::Query options	   
			  'mindate:s'  => \$options{'-mindate'},
			  'maxdate:s'  => \$options{'-maxdate'},
			  'reldate:s'  => \$options{'-reldate'}, 
			  'datetype:s' => \$options{'-datetype'}, # edat or mdat
			  'maxids:i'   => \$options{'-maxids'},
			  'q|query:s'  => \$options{'-query'},
			 );
die "--taxonomy and --primary cannot co-occur" if $taxonomy and $primary;

my $out = Bio::SeqIO->new(-format => $format); # write to STDOUT

my $dbh = Bio::DB::GenBank->new(-verbose => $debug);

sub dump_seq($) {
    my $seq = shift;
    print STDERR "seq->primary_id       = ".$seq->primary_id."\n";
    print STDERR "seq->display_id       = ".$seq->display_id."\n";
    print STDERR "seq->accession_number = ".$seq->accession_number."\n";
    print STDERR "seq->desc             = ".$seq->desc."\n";
}
my %taxonomy_ids;
sub process_id($) {
    $_ = shift;
    print STDERR "\$_ = $_\n" if $debug;
    my ($id, undef) = split;
    print STDERR "\$id = $id\n" if $debug;
    my @id = split /_/, $id;
    die "invalid id for taxonomy: $id" if scalar @id != 2;
    $taxonomy_ids{$id[0]} = $id;
    return $id[0];
}
sub process_seq($) {
    my $seq = shift;
    if ($debug) { print STDERR "process_seq BEFORE\n"; dump_seq($seq); }
    if ($taxonomy) {
        $seq->desc($seq->accession_number." ".$seq->display_id." ".$seq->desc);
        $seq->display_id($taxonomy_ids{$seq->primary_id});
    } elsif ($primary) {
        $seq->desc($seq->accession_number." ".$seq->display_id." ".$seq->desc);
        $seq->display_id($seq->primary_id."_".$seq->display_id);
    }
    print STDERR $seq->display_id."\t".$seq->accession_number."\n" if $accession;
    if ($debug) { print STDERR "process_seq AFTER\n"; dump_seq($seq); }
}
my $query;
if( $gifile ) {
	my @ids;
	open( my $fh => $gifile ) || die $!;
	while(<$fh>) {
        if ($taxonomy) {
		    push @ids, $_;
        } else {
		    push @ids, split;
        }
	}
	close($fh);	
	while( @ids ) {
		my @mini_ids = splice(@ids, 0, $options{'-maxids'});
        @mini_ids = map { process_id($_) } @mini_ids if $taxonomy;
		$query = Bio::DB::Query::GenBank->new(%options,
														  -ids => \@mini_ids,
														 );
		my $stream = $dbh->get_Stream_by_query($query);
		while( my $seq = $stream->next_seq ) {
            process_seq($seq);
			$out->write_seq($seq);
		}
	}
	exit;
} elsif( $options{'-query'}) {
    $options{'-query'} = process_id($options{'-query'}) if $taxonomy;
	$query = Bio::DB::Query::GenBank->new(%options);
} else {
	die("no query string or gifile\n");
}
my $stream = $dbh->get_Stream_by_query($query);
while( my $seq = $stream->next_seq ) {
    process_seq($seq);
	$out->write_seq($seq);
}
