.PHONY: all clean pypi test html readme README.rst inc

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

readme: README.rst

README.rst:
	f=$$(mktemp);\
	  sed '/USAGE-BEGIN/q' README.rst >> $$f; \
	  printf '\n.. code-block::\n' >> $$f;\
	  ./samsift/samsift -h 2>&1 | perl -pe 's/^(.*)$$/\t\1/g' >> $$f; \
	  printf '\n' >> $$f;\
	  sed -n '/USAGE-END/,$$ p' README.rst >> $$f;\
	  cp $$f README.rst

inc:
	./samsift/increment_version.py

test:
	$(MAKE) -C tests
