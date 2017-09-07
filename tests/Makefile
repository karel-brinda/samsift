.PHONY: all clean basic columns id

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

SIFT=../samsift/samsift.py
SV=samtools view

all: basic columns id

basic:
	$(SIFT) -i test.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'
	$(SIFT) -i test.bam -o test_out.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'

id:
	$(SIFT) -i test.bam | grep -v samsift > _test.sam
	diff test.sam _test.sam

columns: test.debug_samtools.sam test.debug_samsift.sam
	diff $^

%.debug_samtools.sam: %.bam
	$(SV) $< \
		| cut -f -11 \
		> $@

%.debug_samsift.sam: %.bam
	($(SIFT) -i $< -d '"\t".join(map(str,[QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]))' \
		> /dev/null) 2>&1 \
		| grep -v samsift \
		| cut -f 3- \
		| cut -f -11 \
		> $@


clean:
	rm -f *.sam test_out.bam
