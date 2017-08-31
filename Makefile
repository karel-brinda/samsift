.PHONY: all clean pypi test

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

all:


clean:
	rm -fr build/ dist/ samsift/__pycache__ samsift.egg-info

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

test:
	$(MAKE) -C test
