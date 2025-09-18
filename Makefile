.PHONY: all clean pypi test test1 test2 readme rst html inc wpypi wconda sha256 yapf

SHELL=/usr/bin/env bash
SAMS=./samsift/samsift.py

.SECONDARY:

all: readme tests/tests/test.bam

clean:
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info
	rm -f *.sam *.bam
	rm -f README.html README.sh
	$(MAKE) -C tests clean

pypi:
	$(MAKE) clean
	$(PYTHON) setup.py sdist bdist_wheel
	$(PYTHON) -m twine upload dist/*

tests/tests/test.bam: tests/tests/test.sam
	samtools view -b "$<" > "$@"

readme: rst html

html:
	rst2html.py README.rst > README.html

rst:
	f=$$(mktemp);\
	  sed '/USAGE-BEGIN/q' README.rst >> $$f; \
	  printf '\n.. code-block::\n' >> $$f;\
	  $(SAMS) -h 2>&1 | perl -pe 's/^(.*)$$/\t\1/g' >> $$f; \
	  printf '\n' >> $$f;\
	  sed -n '/USAGE-END/,$$ p' README.rst >> $$f;\
	  cat $$f \
	  | perl -pe 's/^[\s]+$$/\n/g' \
	  | perl -pe 's/[\s]+$$/\n/g' \
	  > README.rst;

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

wconda:
	open https://bioconda.github.io/recipes/samsift/README.html

wpypi:
	open https://pypi.python.org/pypi/samsift

sha256:
	s=$$(curl https://pypi.python.org/pypi/samsift  2>/dev/null| perl -pe 's/#/\n/g' | grep -o 'https.*\.tar\.gz' | xargs curl -L 2>/dev/null | shasum -a 256 | awk '{print $$1;}'); echo $$s; echo $$s | pbcopy

yapf: ## Run YAPF (inline replacement)
	yapf -i --recursive samsift setup.py tests

