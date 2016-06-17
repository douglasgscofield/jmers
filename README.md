jmers: Kmer-guided processing of jumping reads
==============================================

Clean and process jumping reads using a kmer database built from 'good' sequences of the same genome.  The 'good' sequences are assumed to contain all to nearly all of the kmers found in the target genome, and the jumping reads are assumed to be in various stages of incomplete processing.

Regardless of the protocol used to prepare jumping-read libraries, only the ends of genomic fragments are maintained, with intervening sequence removed and with additional utility sequences added.  Traditional jumping-read processing has involved conservative trim choices due to inability to know with reasonable certainty where the boundaries of the genomic and utility sequences may be found in the fragments ultimately sequences.

Kmer-based processing uses a set of putative 'good' genomic kmers to makes fewer assumptions about what is and is not genomic sequence in jumping-read fragments, and should result in more accurate trimming and, ultimately, longer read lengths after trimming.  Longer lengths of jumping reads enables more accurate mapping and lower link thresholds during scaffolding.

The first target is to process fosmid end sequences.  Some characteristics of the fosmid ends:

1. Production involved a circularisation which ligated the ends of each fosmid together
2. The ligation site is blunt and unknown; it was not created with a restriction enzyme
3. Interior sequence is cut away beginning some distance on either side of the ligation site
4. Utility sequences are attached to these 'inner' ends of the fragment
5. Fosmid-end fragments are sequenced with 250-bp paired-end MiSeq with designed overlap
6. Overlapping with FLASH produces joined fragments around 300-400 bp (will add figure)
7. Conservative trimming resulted in a judgment of keeping the outer 80-90 bp of the end of each joined fragment, after trimming away the outermost 10 bp.
8. We think we can do better than that

**Plan**:

* Read in 'good' genomic sequences (WGS assemblies, superreads, fosmid pool contigs, etc.)
* Create 'good kmer' database
* Read merged fosmid-end fragments
* Walk each fragment, identifying kmers as genomic or non-genomic
* Trim fragment into read 1 and read 2 appropriately
* Output fosmid end read pair
* Get a better genome assembly in the end

**Correcting**: Kmer-walking within known-kmer contexts will be walking through
messy kmers containing sequencing errors.  We might consider doing some simple
correction, perhaps using correction code from some other tool?

**Haplotypes**: Kmer-walking within known-kmer contexts will also be walking
through kmers belonging to haplotypes that are not in the genome assembly code.
Is this distinguishable from correction above?

Both correcting and haplotypes should be addressed after the kmer boundaries
are detected, eller hur?

**Dreaming**:

* Extend to mate-pair processing
* Some form of paired-end contamination removal
* Produce candidate subsets of jumping reads: reads matching 'good kmers' in a sequence subset


Design
------

The namespace is `jmers`.

The `Seq` class (`Seq.h`) is a simple sequence holder with members `name`,
`comment`, `sequence`, `quality`, `l` (sequence length) and `has_quality` and
methods for `fill()`, `write_fastq()`, `write_fasta()`.

The `Input` class (`Input.h`) is a simple sequence reader that uses Heng Li's
[`kseq.h`](http://lh3lh3.users.sourceforge.net/kseq.shtml) to read FastQ- and
Fasta-format input.  The constructor can receive a filename via `std::string`
or `char *`; `read(Seq&)` method fills a `Seq`, returning `false` if
end-of-file; and `close()` shuts things down.  `Input` can currently handle
uncompressed and gzip-compressed input (via `kseq.h`) but will be extended to
handle bzip2 and xz formats and will also be able to read from stdin/fifo.

The `FosmidEndFragment` class (`FosmidEndFragment.h`) initialises with a `Seq`
and implements methods for inferring boundaries, splitting the fragment after
inference into separate `Seq`s for read1 and read2, and writing them out.
Output still needs to be abstracted so the ultimate output format (Fasta,
FastQ, something else) needn't be known by this and similar classes.

Boundary inference might be better abstracted into a separate factory so that
various fragment types needn't duplicate common inference code, for example
sliding contexts from non-genomic to genomic, genomic to chimeric genomic, and
genomic to non-genomic.  That would also make it easier to add additional
context slides from contaminant (kmer db 1) to noncontaminant (kmer db 2), etc.

The `KmerBoundary` class (`KmerBoundary.h`) implements a factory to detect
kmer-determined boundaries between context A and context B.  It initialises
with the type of boundary to detect and the kmer database(s) to use when
detecting the contexts.  It does not initialise with the sense of the boundary,
meaning whether the boundary is detected as a shift from A to B or from B to A.
That is instead specified when the factory is used.

The `KmerBoundary` class is used by passing it a `Seq`, a sense (context A to B
or B to A), and optionally a starting position; the default starting position
is the first position of the `Seq`.

To determine boundaries while handling error-containing reads, we will use a
technique based on that used in
[Quorum](http://www.genome.umd.edu/quorum.html), which is natural since like
Quorum we are using [Jellyfish](http://www.genome.umd.edu/jellyfish.html) to
manage our kmer database.

1. Scan both sides of a boundary against involved kmer databases to determine context to consider.
2. Identify a good 'anchor' kmer to begin accretion
3. Look for an alternative kmer in the desired direction by querying the kmer database for a kmer terminating at each of the four possible bases. We can't choose among these based on kmer depth as quorum does, because we assume our kmers are all good, but we can use backup to retry choices
4. If there is one possibility, then choose it with a minimal penalty
5. If there is more than one possible kmer, accept one and then slide one more base, repeating the process
6. As we slide, if we cross a certain accumulated penalty, then go back and try alternate kmer choices
7. If the accumulated error content exceeds some higher level, reject the kmer and therefore the read
