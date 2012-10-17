TARGETS:=main doc
LDLIBS:=-lgsl -lblas -lstdc++
SRCS:=$(wildcard *.cc)
OBJS:=$(patsubst %cc,%o,$(SRCS))
CXXFLAGS=-g -Wall

default: main

all: depend $(TARGETS)

main: $(OBJS)

.PHONY: doc
doc:
	doxygen doxygen.config

.PHONY: clean
clean:
	@rm -rf $(TARGETS) *.o

depend:
	g++ -MM $(SRCS) > Makefile.dep

include Makefile.dep
# DO NOT DELETE
