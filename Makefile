.PHONY: all clean pypi

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

all:


clean:

pypi:
	/usr/bin/env python3 setup.py sdist bdist_wheel upload

