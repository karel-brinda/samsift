.PHONY: all clean basic columns id seqheadless

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

SIFT=../../samsift/samsift.py
NORM=../../samsift/samsift_norm_sam.py -H
SV=samtools view

all: basic columns id seqheadless

basic:
	$(SIFT) -i test.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'
	$(SIFT) -i test.bam -o _test_out.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'

id:
	$(SIFT) -i test.bam | grep -v samsift > _test.sam
	diff <($(NORM) -i test.sam) <($(NORM) -i _test.sam)

seqheadless:
	$(SIFT) -i seqheadless.bam
columns: _test.debug_samtools.sam _test.debug_samsift.sam
	diff $^

_%.debug_samtools.sam: %.bam
	$(SV) $< \
		| cut -f -11 \
		> $@

_%.debug_samsift.sam: %.bam
	($(SIFT) -i $< -d '"\t".join(map(str,[QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]))' \
		> /dev/null) 2>&1 \
		| grep -v '^\[E::idx_find_and_load\]' \
		| grep -v samsift \
		| cut -f 3- \
		| cut -f -11 \
		> $@


clean:
	rm -f _*
