.PHONY: all clean pypi test test_readme readme README.rst README.html inc

SHELL=/usr/bin/env bash

.SECONDARY:

all: readme


clean:
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info
	rm -f *.sam *.bam

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

readme: README.rst README.html

README.html:
	rst2html.py README.rst > README.html

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
	$(MAKE) readme

test: test_readme
	$(MAKE) -C tests

test_readme:
	cat README.rst \
	     | grep -E '       ' \
	     | perl -pe 's/^\s*//g' \
	     | grep -E "^samsift" \
	     | perl -pe 's/^samsift /samsift\/samsift /g' \
	     | xargs -L 1 bash -x
