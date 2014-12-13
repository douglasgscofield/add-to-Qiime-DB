#!/usr/bin/env perl

# Utility script used for matching up incomplete with complete taxonomic
# hierarchies.  Probably not useful as is, but very useful in conjunction with
# direct inspection of hierarchies.

use strict;
use warnings;
use Getopt::Long;

my ($o_1, $o_2, $o_t, $o_c);
my $verbose = 0;

GetOptions("1=s" => \$o_1,
           "2=s" => \$o_2,
           "taxonomy=s" => \$o_t,
           "completed=s" => \$o_c) or die "unknown option";


my ($n_1, $n_2in1, $n_2, $n_c, $n_cin1, $n_cin2, $n_tin1, $n_tin2) = (0) x 8;
my %T;
# load original taxonomies
my %O1;
my %T1;
open (F1, "<$o_1") or die "cannot open file 1 '$o_1': $!";
while (<F1>) {
    ++$n_1;
    chomp; 
    my @l = split /\t/;
    #unshift @l, "tmp";
    die "ID $l[1] already exists when loading first file" if exists $O1{$l[1]};
    $O1{$l[1]} = $l[2];
    $l[2] =~ / taxonid_(\d+)$/;
    if (not $1) {
        print STDERR "no taxonid found for first file:   $l[1]  $l[2]\n";
    } else {
        $T1{$1} = $l[2];
        $T{$1} = $l[2];
    }
}
print STDERR "Number of hierarchies in first file : $n_1\n";
print STDERR "Number of taxonids in first file    : ".scalar(keys %T1)."\n";
# load larger set of original taxonomies
my %O2;
my %T2;
if ($o_2) {
    open (F2, "<$o_2") or die "cannot open file 2 '$o_2': $!";
    while (<F2>) {
        ++$n_2;
        chomp; 
        my @l = split /\t/;
        #unshift @l, "tmp";
        die "ID $l[1] already exists when loading second file" if exists $O2{$l[1]};
        ++$n_2in1 if exists $O1{$l[1]};
        $O2{$l[1]} = $l[2];
        $l[2] =~ / taxonid_(\d+)$/;
        if (not $1) {
            print STDERR "no taxonid found for second file:   $l[1]  $l[2]\n";
        } else {
            $T2{$1} = $l[2];
            $T{$1} = $l[2];
        }
    }
    print STDERR "Number of hierarchies in second file: $n_2\n";
    print STDERR "Number of taxonids in second file   : ".scalar(keys %T2)."\n";
    print STDERR "Number of these also in first file  : $n_2in1\n";
}

# load completed taxonomies
my %C;
my %TC;
my ($n_just1, $n_just2, $n_both12) = (0) x 3;
my ($n_tjust1, $n_tjust2, $n_tboth12) = (0) x 3;
open (C, "<$o_c") or die "cannot open file of completed taxonomies '$o_c': $!";
while (<C>) {
    ++$n_c;
    chomp; 
    my @l = split /\t/;
    #unshift @l, "tmp";
    die "ID $l[1] already exists when loading completed hierarchies" if exists $C{$l[1]};
    if (exists $O1{$l[1]}) {
        ++$n_cin1;
        ++$n_just1 if not exists $O2{$l[1]};
        #print STDERR "o1: $l[1]    $O1{$l[1]}   $l[2]\n";;
    }
    if (exists $O2{$l[1]}) {
        ++$n_cin2;
        ++$n_just2 if not exists $O1{$l[1]};
        #print STDERR "o2: $l[1]    $O2{$l[1]}   $l[2]\n";;
    }
    ++$n_both12 if exists $O1{$l[1]} and exists $O2{$l[1]};
    $C{$l[1]} = $l[2];
    $l[2] =~ / taxonid_(\d+)$/;
    if (not $1) {
        print STDERR "no-taxonid-found-for-completed\t$l[1]\t$l[2]\t$O1{$l[1]}\t$O2{$l[1]}\tend\n";
    } else {
        if (exists $T1{$1}) {
            ++$n_tin1;
            ++$n_tjust1 if not exists $T2{$1};
            #print STDERR "o1: $l[1]    $O1{$l[1]}   $l[2]\n";;
        }
        if (exists $T2{$1}) {
            ++$n_tin2;
            ++$n_tjust2 if not exists $T1{$1};
            #print STDERR "o2: $l[1]    $O2{$l[1]}   $l[2]\n";;
        }
        if (exists $TC{$1}) {
            if ($verbose) {
                print STDERR "taxonid $1 already exists in \%TC";
                if ($TC{$1} eq $l[2]) {
                    print STDERR "\thierarchies-equivalent";
                } else {
                    print STDERR "\t$TC{$1}\t$l[2]";
                }
                print STDERR "\n";
            }
        }
        if (exists $T1{$1} and exists $T2{$1}) {
            ++$n_tboth12;
            if ($verbose) {
                print STDERR "taxonid $1 in both '$o_1' and '$o_2'";
                if ($T1{$1} eq $T2{$1}) {
                    print STDERR "\thierarchies-equivalent";
                } else {
                    print STDERR "\t$T1{$1}\t$T2{$1}";
                }
                print STDERR "\n";
            }
        }
        $TC{$1} = $l[2];
        print STDOUT "$l[1]\t$1\t$T{$1}\t$TC{$1}\n";
    }
}
print STDERR "Number of completed                 : $n_c\n";
print STDERR "Number of completed in first file   : $n_cin1\n";
print STDERR "Number of completed in second file  : $n_cin2\n";
print STDERR "\n";
print STDERR "Number of accession only first file : $n_just1\n";
print STDERR "Number of accession only second file: $n_just2\n";
print STDERR "Number of accession both files      : $n_both12\n";
print STDERR "\n";
print STDERR "Total accessions                    : ".($n_just1 + $n_just2 + $n_both12)."\n";
print STDERR "\n";
print STDERR "Number of taxonids in first file    : $n_tin1\n";
print STDERR "Number of taxonids in second file   : $n_tin2\n";
print STDERR "\n";
print STDERR "Number of taxonids only first file  : $n_tjust1\n";
print STDERR "Number of taxonids only second file : $n_tjust2\n";
print STDERR "Number of taxonids both files       : $n_tboth12\n";
print STDERR "\n";
print STDERR "Total taxonids                      : ".($n_tjust1 + $n_tjust2 + $n_tboth12)."\n";
print STDERR "\n";
print STDERR "Number of taxonids in files 1 and 2 : ".scalar(keys %T)."\n";
print STDERR "Number of taxonids in completed     : ".scalar(keys %TC)."\n";
print STDERR "\n";

