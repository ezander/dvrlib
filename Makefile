#
TARGETS:=main doc
LDLIBS:=-lgsl -lblas -lstdc++
SRCS:=$(wildcard src/*.cc)
OBJS:=$(patsubst %cc,%o,$(SRCS))
CXXFLAGS=-g -Wall

.PHONY: default
default: depend main

.PHONY: all
all: depend $(TARGETS)

.PHONY: main
main: src/main

src/main: $(OBJS)

.PHONY: doc
doc:
	doxygen doxygen.config

.PHONY: clean
clean:
	@rm -rf $(TARGETS) $(OBJS)

depend:
	@g++ -MM $(SRCS) > Makefile.dep

include Makefile.dep

