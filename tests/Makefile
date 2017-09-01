.PHONY: all clean

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

SIFT=../samsift/samsift

all:
	$(SIFT) -i test.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'
	$(SIFT) -i test.bam -o test_out.bam -c 'ln=len(SEQ);ab=sum(a.query_qualities)/ln'

clean:
	rm -f *.sam
