.PHONY: all help clean cleanall test test_readme test_tests \
        readme rst html \
        format pylint flake8 \
        inc build install sha256 \
        wconda wpypi

SHELL=/usr/bin/env bash -eo pipefail
SAMS=./samsift/samsift.py #used also in tests


PYTHON=/usr/bin/env python3
PIP=/usr/bin/env python3 -m pip

.SECONDARY:

.SUFFIXES:

#################
## BASIC RULES ##
#################

all: readme tests/tests/test.bam

help: ## Print help messages
	@echo -e "$$(grep -hE '^\S*(:.*)?##' $(MAKEFILE_LIST) \
		| sed \
			-e 's/:.*##\s*/:/' \
			-e 's/^\(.*\):\(.*\)/   \\x1b[36m\1\\x1b[m:\2/' \
			-e 's/^\([^#]\)/\1/g' \
			-e 's/: /:/g' \
			-e 's/^#\(.*\)#/\\x1b[90m\1\\x1b[m/' \
		| column -c2 -t -s : )"

clean: ## Clean
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info
	rm -f *.sam *.bam
	rm -f README.html README.sh
	$(MAKE) -C tests clean

cleanall: clean ## Clean all


#############
## TESTING ##
#############


test: ## Test everything
test: test_readme test_tests


tests/tests/test.bam: tests/tests/test.sam
	samtools view -b "$<" > "$@"

test_readme: ## Test README commands
	echo "#! /usr/bin/env bash -e -x" > README.sh
	cat README.rst \
		| grep -E '       ' \
		| perl -pe 's/^\s*//g' \
		| grep -E "^samsift" \
		| perl -pe 's@^samsift/samsift@samsift@g' \
		| perl -pe 's@^samsift @$(SAMS) @g' \
		>> README.sh
	/usr/bin/env bash README.sh

test_tests: ## Run all tests
	$(MAKE) -C tests


#################
## MAINTENANCE ##
#A###############

format: ## Run YAPF (inline replacement)
	yapf -i --recursive samsift .py tests

pylint: ## Run PyLint
	$(PYTHON) -m pylint attotree

flake8: ## Run Flake8
	flake8


###############
## PACKAGING ##
###############

inc: ## Increment version
	./samsift/increment_version.py
	$(MAKE) readme

build: ## Build
	$(PYTHON) -m build
	$(PYTHON) -m twine check dist/*


install: ## Install using PIP (current env)
	$(PIP) uninstall -y samsift || true
	$(PIP) install .


sha256: ## Compute sha256 for the PyPI package
	s=$$(curl https://pypi.python.org/pypi/samsift  2>/dev/null| perl -pe 's/#/\n/g' | grep -o 'https.*\.tar\.gz' | xargs curl -L 2>/dev/null | shasum -a 256 | awk '{print $$1;}'); echo $$s; echo $$s | pbcopy


#########################
## DOCUMENTATION & WEB ##
#########################

readme: ## Update README files
readme: rst html

html: ## Convert README to HTML
	rst2html.py README.rst > README.html

rst: ## Update help message in README
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

wconda: ## Open Bioconda webpage
	open https://bioconda.github.io/recipes/samsift/README.html

wpypi: ## Open PyPI webpage
	open https://pypi.python.org/pypi/samsift

