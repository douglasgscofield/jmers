#!/usr/bin/env perl

use strict;
use warnings;
use feature 'say';
use Getopt::Long;
use List::Util;
use Math::Random;

# NOTES:
# joinpos at which reads are butt-joined suffixes the read name as _j:pos and is 0-based
#
# jellyfish database from a single fosmid sequence
# from this same sequence
# simulate 250-bp read ends definitely not overlapping
# butt-join pairs into 500-bp fragments
# randomly trim fixed total of 100bp from both ends of fragment, tracking join location
#   pick length of left end from N(mean=50,sd=10) and set right end to 100 - left
#   keep track of join location
#
# when later adding errors, something like:
# ./wgsim -N10 -s0 -e0 -r0 -d400 -1250 -2250 -R0 fosmids/999.fasta r1.fq r2.fq
# wgsim overlapping reads from each segment

my $ref;
my ($read1, $read2) = qw/ r1.fq r2.fq /;
my $o_N = 5;  # number of read pairs to simulate
my $o_pair_outer = 1000;
my $o_pair_se = 0;
my $o_readlen = 250;
my $o_trim = 100;
my $o_trim_sd = 20;
my $o_seed = 1;
GetOptions(
    "N=i" => \$o_N,
) or die "unknown option: $!";
$ref = shift if @ARGV;
die "cannot locate reference fasta $ref" if ! -f $ref;

# Simulate non-overlapping read ends
my $cmd = "./wgsim/wgsim -S$o_seed -N$o_N -d$o_pair_outer -s$o_pair_se -1$o_readlen -2$o_readlen -e0 -R0 $ref $read1 $read2 > /dev/null";
say STDERR "wgsim command: $cmd";
system($cmd);

# Merge together
open(my $r1, "<", $read1) or die "could not open $read1: $!";
open(my $r2, "<", $read2) or die "could not open $read2: $!";

while (<$r1>) {
    my $r1_1 = $_;    chomp $r1_1; die "unknown readname format $r1_1" if $r1_1 !~ /\/1$/;
    my $r1_2 = <$r1>; chomp $r1_2;
    my $r1_3 = <$r1>; chomp $r1_3;
    my $r1_4 = <$r1>; chomp $r1_4;
    my $r2_1 = <$r2>; chomp $r2_1; die "unknown readname format $r2_1" if $r2_1 !~ /\/2$/;
    my $r2_2 = <$r2>; chomp $r2_2;
    my $r2_3 = <$r2>; chomp $r2_3;
    my $r2_4 = <$r2>; chomp $r2_4;

    my $s = $r1_2 . $r2_2;
    my $q = $r1_4 . $r2_4;
    my $joinpos = $o_readlen; # note this is 0-based
    my $lefttrim = sprintf "%.0f", Math::Random::random_normal(1, $o_trim/2, $o_trim_sd);
    my $righttrim = $o_trim - $lefttrim;
    die "some trim length negative" if $lefttrim < 0 or $righttrim < 0;
    $s = substr $s, $lefttrim; $s = substr $s, 0, -$righttrim;
    $q = substr $q, $lefttrim; $q = substr $q, 0, -$righttrim;
    $joinpos -= $lefttrim;
    #$joinpos += 1; # if uncommented, converts to 1-based

    my $m_1 = $r1_1; $m_1 =~ s/\/1$/_j:$joinpos/;
    my $m_2 = $s;
    my $m_3 = "+";
    my $m_4 = $q;

    print STDOUT "$m_1\n$m_2\n$m_3\n$m_4\n";
}
