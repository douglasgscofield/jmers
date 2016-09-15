JELLYFISH= ../../Jellyfish/bin/jellyfish
.SUFFIXES:.fasta .fa  .jf


.fasta.jf:
	$(JELLYFISH) count -m 21 -s 10G -t 1 $< --out-counter-len 2 -o $@

.fa.jf:
	$(JELLYFISH) count -m 21 -s 10G -t 1 $< --out-counter-len 2 -o $@

# initiates the DWGSIM submodule
.init:
	cd DWGSIM; git submodule init; git submodule update; make
