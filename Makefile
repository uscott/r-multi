RCMD := R CMD

.DEFAULT_GOAL := all

all: clean check build
.PHONY: all

build:
	$(RCMD) build .
.PHONY: build

check:
	$(RCMD) check .
.PHONY: check

clean:
	rm -f *.tar.gz
	rm -f src/*.o src/*.so
.PHONY: clean

format:
	clang-format -i src/*.c src/*.h
.PHONY: format
