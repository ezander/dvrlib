#
TARGETS:=main doc
LDLIBS:=-lgsl -lstdc++ -Lgsl/.libs\
	-Lgsl/cblas/.libs -lgslcblas
SRCS:=$(wildcard src/*.cc)
OBJS:=$(patsubst %cc,%o,$(SRCS))
CXXFLAGS=-g -Wall -Igsl -Igsl/cblas

.PHONY: default
default: depend main lib


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
	@rm -rf $(TARGETS) $(OBJS) libdvrlib.a

depend:
	@g++ -MM $(SRCS) > Makefile.dep

.PHONY: lib
lib:
	ar rsv libdvrlib.a src/*.o

include Makefile.dep

