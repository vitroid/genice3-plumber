.PHONY: help test

help:
	@echo "Targets:"
	@echo "  make test    Run smoke tests (test/test1.sh)"

test:
	poetry run bash test/test1.sh
