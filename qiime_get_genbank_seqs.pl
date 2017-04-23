#!/usr/bin/env perl 

=head1 NAME

qiime_get_genbank_seqs.pl - Retrieve Genbank sequences and produce output useful
                            for augmenting a Qiime reference database.
 
Two files are output: a file of FASTA-format sequences to STDOUT, and a file of
IDs to STDERR.

=head1 USAGE

qiime_get_genbank_seqs.pl --accessionfile seqs.ids > seqs.fa 2>seqs_ids.txt

This script is derived from the BioPerl script bp_download_query_genbank.pl, and
is modified to support its use with augmenting a Qiime reference database as
part of the scripts in L<< add-to-Qiime-DB|https://github.com/douglasgscofield/add-to-Qiime-DB >>.

Queries in the accession file are in the form of ACCESSION.VERSION <tab> TAXONID.
Each downloaded sequence is given the name ACCESSION.VERSION_taxTAXONID,
with '_tax' as a separator between the versioned accession number and the taxon ID.

=head2 Query options

These options modify the manner in which queries are submitted to GenBank.

=over

=item --maxids maximum number of IDs to retrieve in a set (default 100)

=item --reldate 

=item --maxdate maxdate for a record

=item --mindate minimum date for record

=item --datetype edat or mdat (entered or modified)

=back

=head1 AUTHOR Jason Stajich

Jason Stajich, jason-AT-bioperl.org

Heavily modified by Douglas Scofield (Uppsala University, douglas.scofield@ebc.uu.se) for Qiime-oriented work as part of

    L<< add-to-Qiime-DB|https://github.com/douglasgscofield/add-to-Qiime-DB >>

=cut

use strict;
use warnings;
use Bio::DB::GenBank;
use Bio::DB::GenPept;
use Bio::DB::Query::GenBank;
use Bio::SeqIO;
use Getopt::Long;
use Pod::Usage;

my ( $format, $debug, $primary, $taxsep, %options );

# for Qiime work, all four of these should remain unchanged
$format         = 'fasta';
$taxsep         = '_tax';
$options{'-db'} = 'nucleotide';

$options{'-maxids'} = '100';
my $accessionfile;
GetOptions(
    'h|help' => sub {
        exec( 'perldoc', $0 );
        exit(0);
    },
    'v|verbose'       => \$debug,
    'accessionfile:s' => \$accessionfile,

    # DB::Query options
    'mindate:s'  => \$options{'-mindate'},
    'maxdate:s'  => \$options{'-maxdate'},
    'reldate:s'  => \$options{'-reldate'},
    'datetype:s' => \$options{'-datetype'},    # edat or mdat
    'maxids:i'   => \$options{'-maxids'},
    'q|query:s'  => \$options{'-query'},
);

my $out = Bio::SeqIO->new( -format => $format );    # write to STDOUT

my $dbh = Bio::DB::GenBank->new( -verbose => $debug );

sub dump_seq($) {
    my $seq = shift;
    print STDERR "seq->primary_id       = " . $seq->primary_id . "\n";
    print STDERR "seq->display_id       = " . $seq->display_id . "\n";
    print STDERR "seq->accession_number = " . $seq->accession_number . "\n";
    print STDERR "seq->version          = " . $seq->version . "\n";
    print STDERR "seq->desc             = " . $seq->desc . "\n";
}
my %taxonomy_ids;

sub process_id($) {
    $_ = shift;
    print STDERR "\$_ = $_\n" if $debug;
    my ( $acc, $taxid, undef ) = split;
    print STDERR "\$acc = $acc, \$taxid = $taxid\n" if $debug;
    $taxonomy_ids{ $acc } = $acc . $taxsep . $taxid;
    return $acc;
}

sub process_seq($) {
    my $seq = shift;
    if ($debug) { print STDERR "process_seq BEFORE\n"; dump_seq($seq); }
    my $acc_v = $seq->accession_number . "." . $seq->version;
    $seq->desc( $acc_v . " " . $seq->display_id . " " . $seq->desc );
    $seq->display_id( $taxonomy_ids{ $acc_v } );
    print STDERR $seq->display_id . "\t" . $acc_v . "\n";
    if ($debug) { print STDERR "process_seq AFTER\n"; dump_seq($seq); }
}
my $query;
if ($accessionfile) {
    my @ids;
    open( my $fh => $accessionfile ) || die $!;
    while (<$fh>) {
        push @ids, $_;
    }
    close($fh);
    while (@ids) {
        my @mini_ids = splice( @ids, 0, $options{'-maxids'} );
        @mini_ids = map { process_id($_) } @mini_ids;
        $query = Bio::DB::Query::GenBank->new( %options, -ids => \@mini_ids, );
        my $stream = $dbh->get_Stream_by_query($query);
        while ( my $seq = $stream->next_seq ) {
            process_seq($seq);
            $out->write_seq($seq);
        }
    }
    exit;
}
else {
    die("no accessionfile\n");
}
my $stream = $dbh->get_Stream_by_query($query);
while ( my $seq = $stream->next_seq ) {
    process_seq($seq);
    $out->write_seq($seq);
}
