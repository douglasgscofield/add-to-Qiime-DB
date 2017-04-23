#!/usr/bin/env perl

# Parse blast results for GenBank IDs and hit coordinates.  Part of
#
#      https://github.com/douglasgscofield/add-to-Qiime-DB
#
# Copyright (c) 2016, Douglas G. Scofield
# Evolutionary Biology Centre, Uppsala University
# douglas.scofield@ebc.uu.se

use strict;
use warnings;
use Getopt::Long;

my $o_with_gi = 0;
my $o_help = 0;

sub usage() {
print STDERR "
USAGE

    $0 [ --with-gi ] blast-output-table > id-table

If the blast table uses subject IDs containing GIs (prior to Blast 2.5.0+),
specify the --with-gi option.

Use -h or --help for this usage message.

OUTPUT TO STDOUT

A table containing the versioned accession number and the taxon ID, separated
by an underscore, followed by the start and end of the hit in the target.

    accession.version_taxon-ID <tab> start <tab> end

With --with-gi, the subject ID is expected to be in the previous NCBI format,
which contains 'gi|GI-number|gb|accession.version'.  A somewhat different
format is written to stdout, with the GI number and taxon IDs separated by an
underscore, and the version is stripped off the accession number.

    GI-number_taxon-ID <tab> accession <tab> start <tab> end

";
    exit(1);
}

GetOptions("with-gi" => \$o_with_gi, "h|help" =>\$o_help) and !$o_help or usage();

while (<>) {
    chomp;
    my @l = split /\t/;
    if ($o_with_gi) {
        my @ll = split /\|/, $l[1];    # parse the target ID returned by Blast gi|gi-number|gb|accession.version|
        my $gi = $ll[1];               # the gi ID
        my $accession = $ll[3];        # accession.version
        $accession =~ s/\.[0-9]+$//;   # strip version number off the end: .1, .2, etc.
        my $start = $l[8];             # start of hit in target: > end if - strand
        my $end   = $l[9];             # start of hit in target: < start if - strand
        my $tid   = $l[12];            # taxonid
        print STDOUT "${gi}_${tid}\t$accession\t$start\t$end\n";
    } else {
        my $accession = $l[1];         # accession.version
        my $start = $l[8];             # start of hit in target: > end if - strand
        my $end   = $l[9];             # start of hit in target: < start if - strand
        my $tid   = $l[12];            # taxonid
        print STDOUT "${accession}_${tid}\t$start\t$end\n";
    }
}

