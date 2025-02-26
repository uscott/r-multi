RCMD := R CMD
RM_OBJ_FILES := rm -f pkg/src/*.o pkg/src/*.so

.DEFAULT_GOAL := all

all: clean check build
.PHONY: all

build:
	$(RCMD) build pkg
.PHONY: build

check:
	$(RCMD) check pkg
	$(RM_OBJ_FILES)
.PHONY: check

clean:
	rm -f *.tar.gz
	$(RM_OBJ_FILES)
.PHONY: clean

format:
	clang-format -i pkg/src/*.c pkg/src/*.h
.PHONY: format

setup:
	pre-commit install
.PHONY: setup
