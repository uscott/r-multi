RCMD := R CMD

.DEFAULT_GOAL := all

all: clean build
.PHONY: all

build:
	$(RCMD) build .
.PHONY: build

check:
	$(RCMD) check .
.PHONY: check

clean:
	rm -f *.tar.gz
.PHONY: clean
