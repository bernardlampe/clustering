CC := clang++
CCFLAGS := -g -O2 -Wall -Wextra -I.. -std=c++11

SRCS := $(wildcard *.cc)
BINS := $(basename $(SRCS))

.PHONY: all
all: $(BINS)

%: %.cc
	$(CC) $(CCFLAGS) $< -o $@

.PHONY: test
test: all
	@rc=0; \
	for b in $(BINS); do \
		echo "== $$b =="; \
		./$$b || rc=1; \
	done; \
	exit $$rc

.PHONY: clean
clean:
	@rm -rf $(BINS) $(addsuffix .dSYM, $(BINS))
