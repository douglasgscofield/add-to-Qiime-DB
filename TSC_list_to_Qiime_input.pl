#!/usr/bin/env perl

# Utility script to take output files from TSC ITS clustering and produce a file
# suitable for input to Qiime for creation of consensus sequences.  Part of
#
#      https://github.com/douglasgscofield/add-to-Qiime-DB
#
# For more information and an example of use, see the repository README.md.
#
# Copyright (c) 2016, Douglas G. Scofield
# Evolutionary Biology Centre, Uppsala University
# douglas.scofield@ebc.uu.se

use strict;
use warnings;
use Getopt::Long;

my $fasta;
my $otu_list;

GetOptions("fasta=s" => \$fasta,
           "otulist=s" => \$otu_list)
or die "unknown option";
die "must supply both --fasta and --otulist options" if not $fasta or not $otu_list;

open (F, "<$fasta") or die "could not open $fasta: $!";
open (O, "<$otu_list") or die "could not open $otu_list $!";

my %F;
my $F_index = 0;
while (<F>) {
    if (/^>/) {
        chomp;
        my $name = substr($_, 1);
        $F{$F_index++} = $name;
    }
}
die "key and index inconsistency" if scalar(keys %F) != $F_index;
close(F);

while (<O>) {
    chomp;
    my @l = split /\t/;
    my @s = split /,/, $l[2];
    my @n = map { $F{$_} } @s;
    print STDOUT $l[0], "\t", join("\t", @n), "\n";
}
