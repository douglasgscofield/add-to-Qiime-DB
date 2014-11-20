#!/usr/bin/env perl

use strict;
use warnings;

while (<>) {
    chomp;
    my @l = split;
    my $tid = $l[12];
    my @ll = split /\|/, $l[1];
    my $gi = $ll[1];
    my $gb = $ll[3];
    $gb =~ s/\.[0-9]+$//;
    #my $tid = $l[1];
    print STDOUT "${gi}_${tid}\t$gb\n";
}
