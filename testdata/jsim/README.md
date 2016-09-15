Simulating read data
====================

To test jmers' algorithms, simulate paired-end read data for a known reference
sequence with known error rates.  Join the reads into a single read with a
known, randomised join position.  Then, create a Jellyfish kmer database for
the reference and attempt to identify the join position within the read.

I've made Nils Homer's `dwgsim` a submodule here.  To build:

    make .init

A Jellyfish database is not used here, but an appropriate one can be built for
the reference sequence in use (its filename must end with `.fasta` or `.fa`) with

    make ref.jf

For a reference, I've chosen the fosmid scaffold 999.fasta from the Cossu et
al. dataset S4, containing a single 41853 bp sequence with no N.  The `999.jf`
Jellyfish database is for this sequence.  See below for more.


Error-free reads
----------------

As `jsim.pl` is configured now, it simulates the fragment assuming error-free
non-overlapping read pairs.  Read pairs are simulated, and then butt-joined
into a single fragment.  Then, *t*=100 bases are trimmed from the ends of the
fragment: a normal deviate *d* of mean *t*/2 and standard deviation 10 is drawn
and rounded and this many bases are trimmed from the left end of the fragment,
with *t* - *d* bases trimmed from the right end. The trim amount is truncated
of it would exceed [0, *t*]. The resulting fragment is 400 bp in length.

The position of the butt-join is tracked through this trim process and the
fragment is output using a read name having the `/1` of the first readname
replaced with `_j:###` where `###` is the position within the sequence
(0-based) where the second read starts.

To simulate 100 reads, error-free:

    ./jsim.pl -N 100 999.fasta > 999_joined_e0.fq

Two example simulated reads:

```
@999_27610_28609_0:0:0_0:0:0_7_j:197
TGTTTTTTCCCCCTTGGGGGTTTTGATTCCCTTGGTTCTGGTTGCCCCCAGTAGGTTGGTTCTGTTGGGGGTGATTCCCGCCTTGACCTGGGTTGCTATTGTTCCAGTTTTTTCCTTTATTTTTGGCATTATTCGATTGATGGGTTAGATTTTTCATATTCAACTGTTCAAAATTTTGATATTGCTAAATGGGAGAGGGCCAAATCATTTTTGACCTTCTTCCAATTACCAATCCACCATGACAATGGTCTTGAACTACTTTCTAATTTTAAACAAACCTTAGCCACACATATCACTGATCACATTCACGAGTGGCGTCGCCGACGTAGTTTGTGAAAAGCAGAAACCACCAAACAATAATGTCTTGATTGGTTTCTCAAATCACTTTTCTCTCTCCTCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@999_10527_11526_0:1:0_0:0:0_8_j:189
AATGTAGCTGAATTCTAGCTTCAATTATTTTCCAACATCCTCGAGAATGATCAATTAATCGAGAAGTTGCACAATGCGTTATACCTAGTTTCAGGGAATGAAGAACCCATAGAAGCAGGAATAAGGTCTAATCTTCCTGAGAAGAGTACTCAGAGAGTGCAGAGAAGAAAGGTACGAAATGGAAAATAAATCCATCTGCATCTCTACTTCCAACAGTTCCTTGTTTTGAAACCTATACTTGATGGCTTGCTTCACATCGTCCCAATCTCTGATTAAGGCTTTATGATTGGCCCACCATCTGACAAGGGTGTCCTGAAAGGCTACATTTAGTACTGAGATCTTTTGATCCTCTACTACCTTTTCGTATGAGATTTGTTGCCAGCTATCAAATGCCAGTTCA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

```


Reads with base substitution errors
-----------------------------------

Specify a `-e` value.  We can simulate mutations in the reference sequence with
`-r`, but that is not (yet) what we want for `jmers` testing.

    ./jsim.pl -e 0.01 999.fasta > 999_joined_e0.01.fq


Building a Jellyfish database
-----------------------------

    make 999.jf

Using the following gmake Makefile.  Should the database be canonical?  Any other options?

    JELLYFISH= ../../Jellyfish/bin/jellyfish
    .SUFFIXES:.fasta .fa  .jf

    .fasta.jf:
        $(JELLYFISH) count -m 21 -t 1 $< --out-counter-len 2 -o $@

    .fa.jf:
        $(JELLYFISH) count -m 21 -t 1 $< --out-counter-len 2 -o $@


Fosmid end orientation?
-----------------------

We can sort out the orientation of sequences on the joined fosmid ends, but
what is the original orientation of the ends with respect to the genome?  For
a joined fosmid end fragment with an identified split between `AB`:

    A>A>A>AB>B>B>B

which is the original genome orientation and complement orientation (lowercase
is reverse-complement)?

    A>A>A>A   B>B>B>B
    A>A>A>A   b<b<b<b
    a<a<a<a   b<b<b<b
    a<a<a<a   B>B>B>B

The circularisation implies an orientation, of course, but I can't sort it out
just now.  The default is 'mate pair' orientation, which is the first of the
four.  With `-o 0`, the second is used.

