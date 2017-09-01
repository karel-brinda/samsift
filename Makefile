.PHONY: all clean pypi test html desc inc

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

desc:
	./samsift/samsift -h 2>&1 | pbcopy

inc:
	./samsift/increment_version.py

test:
	$(MAKE) -C tests
