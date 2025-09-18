.PHONY: all clean basic columns id

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

#don't keep warnings (they go to the saved output, breaking diffs, and different PySAM versions use different warnings)
SIFT=HTS_LOG_LEVEL=critical ../../samsift/samsift.py
NORM=../../samsift/samsift_norm_sam.py -H
SV=samtools view

all: basic columns id

basic:
	$(SIFT) -i test.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'
	$(SIFT) -i test.bam -o _test_out.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'

id:
	$(SIFT) -i test.bam | grep -v samsift > _test.sam
	diff <($(NORM) -i test.sam) <($(NORM) -i _test.sam)

columns: _test.debug_samtools.sam _test.debug_samsift.sam
	diff $^

_%.debug_samtools.sam: %.bam
	$(SV) $< \
		| cut -f -11 \
		> $@

_%.debug_samsift.sam: %.bam
	($(SIFT) -i $< -d '"\t".join(map(str,[QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL]))' \
		> /dev/null) 2>&1 \
		| grep -v samsift \
		| cut -f 3- \
		| cut -f -11 \
		> $@


clean:
	rm -f _*
