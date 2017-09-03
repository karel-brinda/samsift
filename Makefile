.PHONY: all clean pypi test test1 test2 readme README.rst README.html inc

SHELL=/usr/bin/env bash
SAMS=./samsift/samsift.py

.SECONDARY:

all: readme


clean:
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info
	rm -f *.sam *.bam
	rm -f README.html README.sh
	$(MAKE) -C tests clean

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

readme: README.rst README.html

README.html:
	rst2html.py README.rst > README.html

README.rst:
	f=$$(mktemp);\
	  sed '/USAGE-BEGIN/q' README.rst >> $$f; \
	  printf '\n.. code-block::\n' >> $$f;\
	  $(SAMS) -h 2>&1 | perl -pe 's/^(.*)$$/\t\1/g' >> $$f; \
	  printf '\n' >> $$f;\
	  sed -n '/USAGE-END/,$$ p' README.rst >> $$f;\
	  cat $$f \
	  | perl -pe 's/^[\s]+$$/\n/g' > README.rst;

inc:
	./samsift/increment_version.py
	$(MAKE) readme

test: test1 test2

test1:
	echo "#! /usr/bin/env bash -e -x" > README.sh
	cat README.rst \
		| grep -E '       ' \
		| perl -pe 's/^\s*//g' \
		| grep -E "^samsift" \
		| perl -pe 's@^samsift/samsift@samsift@g' \
		| perl -pe 's@^samsift @$(SAMS) @g' \
		>> README.sh
	bash README.sh

test2:
	$(MAKE) -C tests
