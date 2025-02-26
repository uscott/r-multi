RCMD := R CMD

.DEFAULT_GOAL := all

all: clean check build
.PHONY: all

build:
	$(RCMD) build pkg
.PHONY: build

check:
	$(RCMD) check pkg
.PHONY: check

clean:
	rm -f *.tar.gz
	rm -f pkg/src/*.o pkg/src/*.so
.PHONY: clean

format:
	clang-format -i pkg/src/*.c pkg/src/*.h
.PHONY: format
