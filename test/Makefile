.PHONY: all clean

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

SIFT=../samsift/samsift.py

all:
	$(SIFT) -i test.bam

clean:
	rm -f *.sam
