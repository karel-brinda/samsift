.PHONY: all clean pypi test html

SHELL=/usr/bin/env bash

.SECONDARY:

all: html


clean:
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info
	rm -f *.sam *.bam

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

html:
	rst2html.py README.rst > README.html

test:
	$(MAKE) -C tests
