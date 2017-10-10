.PHONY: all clean filt

SF=../../samsift/samsift.py

SHELL=/usr/bin/env bash -euc -o pipefail
BAM=_test.bam
SMALLBAM=_test.small.bam

.SECONDARY:

all: filt

filt: $(BAM)
	/usr/bin/env python3 -m cProfile -o _filt.prof $(SF) -i "$<" -f 'NM>0' > /dev/null

$(BAM):
	wget -O $(BAM) ftp://ftp.sra.ebi.ac.uk/vol1/ERZ004/ERZ004000/accepted_hits.bam

$(SMALLBAM): $(BAM)
	samtools view -h $< | head -n 10000 | samtools view -b > $@

clean:
	rm -f _*

