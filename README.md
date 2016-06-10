jmers: Kmer-guided processing of jumping reads
==============================================

Clean and process jumping reads using a kmer database built from 'good' sequences of the same genome.  The 'good' sequences are assumed to contain all to nearly all of the kmers found in the target genome, and the jumping reads are assumed to be in various stages of incomplete processing.

Regarldess of the protocol used to prepare jumping-read libraries, only the ends of genomic fragments are maintained, with intervening sequence removed and with additional utility sequences added.  Traditional jumping-read processing has involved conservative trim choices due to inability to know with reasonable certainty where the boundaries of the genomic and utility sequences may be found in the fragments ultimately sequences.

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

Plan:

* Read in 'good' genomic sequences (WGS assemblies, superreads, fosmid pool contigs, etc.)
* Create 'good kmer' database
* Read merged fosmid-end fragments
* Walk each fragment, identifying kmers as genomic or non-genomic
* Trim fragment into read 1 and read 2 appropriately
* Output fosmid end read pair
* Get a better genome assembly in the end


![Plan](plan_20160610.jpg)


Dreaming:

* Extend to mate-pair processing
* Some form of paired-end contamination removal
* Produce candidate subsets of jumping reads: matching 'good kmers' in a sequence subset
