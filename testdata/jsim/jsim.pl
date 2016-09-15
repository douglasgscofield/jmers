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
# TODO?: simulate overlapping MiSeq reads and use FLASH

my $o_N = 1;  # number of read pairs to simulate
my $o_pair_outer = 1000;
my $o_pair_sd = 0;
my $o_readlen = 250;
my $o_trim = 100;
my $o_trim_sd = 10;
my $o_base_error_rate = 0;
my $o_mutation_rate = 0;
my $o_seed = 1;
my $o_orientation = 1;
my $o_wgsim = 0;

my $N_trunc = 0;
my $usage = "
USAGE:  $0  [ options ] ref.fasta > simreads_joined.fq

Simulated joined read fragments are written to stdout in FastQ format.  Under
most circumstances we want to specify sequencing error rate with -e, not
reference mutation rate with -r.

  -N INT     number of read pairs/fragments to simulate [default $o_N]
  -e FLOAT   per-base sequencing error rate [default $o_base_error_rate]
  -r FLOAT   per-base mutation to apply to reference [default $o_mutation_rate]
  -d INT     read pair outer insert size [default $o_pair_outer]
  -s INT     read pair insert size standard deviation [default $o_pair_sd]
  -l INT     read length (applies to both reads of a pair) [default $o_readlen]
  -t INT     total bases to trim from both ends of read fragment [default $o_trim]
  -T INT     half-trim standard deviation [default $o_trim_sd]
  -o 0 or 1  read orientation 0: F,R, 1: F,F [dwgsim only, default $o_orientation]
  -w         use wgsim (https://github.com/lh3/wgsim) to simulate reads; default
             is dwgsim (https://github.com/nh13/DWGSIM).  Both are submodules.

The end trim (controlled with -t and -T) trims t bases in total from the ends
of the fragment.  The left trim is t_L=Normal(mean=t/2,sd=T) and the right trim
is t_R=t-t_L.  The net effect is to produce consistently sized fragments for
which the join location is randomised.  If t_L<0 or t_L>t, it is truncated; a
message is produced if more than 1% of trims are truncated.

No indels are simulated (-R0 to dwgsim and wgsim).

";

GetOptions(
    "d|pair-outer=i"        => \$o_pair_outer,
    "e|base-error-rate=f"   => \$o_base_error_rate,
    "l|read-len=i"          => \$o_readlen,
    "N=i"                   => \$o_N,
    "o|orientation=i"       => \$o_orientation,
    "r|mutation-rate=f"     => \$o_mutation_rate,
    "s|pair-sd=i"           => \$o_pair_sd,
    "t|trim=i"              => \$o_trim,
    "T|trim-sd=i"           => \$o_trim_sd,
    "w|wgsim"               => \$o_wgsim,
) or die "unknown option: $!";
die "$usage" if ! @ARGV;
my $ref = shift @ARGV;
die "cannot locate reference fasta $ref" if ! -f $ref;
my ($prefix, $read1, $read2);
if ($o_wgsim) {
    ($read1, $read2) = ("read_wgsim.read1.fastq",      "read_wgsim.read2.fastq");
} else {
    $prefix = "read_dwgsim";
    ($read1, $read2) = ("$prefix.bwa.read1.fastq", "$prefix.bwa.read2.fastq");
}
#
# Simulate non-overlapping read ends
my $cmd;
if ($o_wgsim) {
    $cmd = "./wgsim/wgsim -S$o_seed -N$o_N -d$o_pair_outer -s$o_pair_sd -1$o_readlen -2$o_readlen -e$o_base_error_rate -R0 -r$o_mutation_rate $ref $read1 $read2 > /dev/null";
} else {
    $cmd = "./DWGSIM/dwgsim -z$o_seed -N$o_N -d$o_pair_outer -s$o_pair_sd -1$o_readlen -2$o_readlen -e$o_base_error_rate -E$o_base_error_rate -R0 -r$o_mutation_rate -y0 -S$o_orientation $ref $prefix > /dev/null";
}
say STDERR "sim command: $cmd";
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
    if    ($lefttrim < 0)       { $lefttrim = 0;       ++$N_trunc; }
    elsif ($lefttrim > $o_trim) { $lefttrim = $o_trim; ++$N_trunc; }
    my $righttrim = $o_trim - $lefttrim;
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

say STDERR "*** $N_trunc of $o_N trim lengths were truncated to keep trim to $o_trim bases" if ($N_trunc / $o_N) > 0.01;

