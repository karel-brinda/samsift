.PHONY: all

SHELL=/usr/bin/env bash -euc -o pipefail

.SECONDARY:

all:
	$(MAKE) -C tests

%:
	$(MAKE) -C tests $@

