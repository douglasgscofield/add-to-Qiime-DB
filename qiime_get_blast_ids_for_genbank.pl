#!/usr/bin/env perl

# Parse blast results for GenBank IDs and hit coordinates.  Part of
#
#      https://github.com/douglasgscofield/add-to-Qiime-DB
#
# Copyright (c) 2016, 2017 Douglas G. Scofield
# Evolutionary Biology Centre, Uppsala University
# douglas.scofield@ebc.uu.se

use strict;
use warnings;
use Getopt::Long;

my $o_help = 0;

sub usage() {
print STDERR "
USAGE

    $0  blast-output-table  >  id-table

If the subject IDs in the blast table contain GIs, they will be ignored and the
versioned accession number will be used instead.  NCBI is phasing out the use
of GIs.

Output to STDOUT is a table containing the versioned accession number, the
taxon ID, and the start and end of the hit in the target.

    accession.version <tab> taxon-ID <tab> start <tab> end

";
    exit(1);
}

GetOptions("h|help" =>\$o_help) and !$o_help or usage();

my $gi_seen = 0;

while (<>) {
    chomp;
    my @l = split /\t/;
    my $accession_version;
    if ($l[1] =~ /^gi\|/) {
        ++$gi_seen;
        # parse the target ID returned by Blast gi|gi-number|gb|accession.version|
        (undef, undef, undef, $accession_version, undef) = split /\|/, $l[1];
    } else {
        $accession_version = $l[1];
    }
    my ($taxonid, $start, $end) = ($l[12], $l[8], $l[9]);
    print STDOUT "$accession_version\t$taxonid\t$start\t$end\n";
}

