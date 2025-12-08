CC := clang++
CCFLAGS := -g -Wall -I.. -std=c++11

SRCS := $(wildcard *.cc)
BINS := $(basename $(SRCS))

.PHONY: all
all: $(BINS)

%: %.cc
	$(CC) $(CCFLAGS) $< -o $@

.PHONY: clean
clean:
	@rm -rf $(BINS) $(addsuffix .dSYM, $(BINS))

