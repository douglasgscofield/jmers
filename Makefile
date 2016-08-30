#
# Basic Makefile for jmers. Currently it doesn't work on OS X, but that will be 
# a later implementation.
#
# This file is heavily based on Jellyfish/examples/query_per_sequence/Makefile.
#
# By: Martin Norling, 2016
#

CC = g++
CXXFLAGS = -Wall -O3
LDFLAGS = -lz 

# build yaggo, add it to PATH so it can be found during
YAGGO_BASE = yaggo
YAGGO = $(YAGGO_BASE)/yaggo
# build Jellyfish, prefix is build/, make and make install
# autoreconf -i
# ./configure --prefix=build/
JELLYFISH_BASE = Jellyfish/build

ifeq ($(OS),Windows_NT)
    CC = g++
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Darwin)
        #CC = clang++
        CC = g++
        #
        # Default values for OS X, since there's usually no pkg-config
        #
        CXXFLAGS += -I$(JELLYFISH_BASE)/include/jellyfish-2.2.6 -std=c++11
        LDFLAGS += -L$(JELLYFISH_BASE)/lib/ -ljellyfish-2.0
    else
        CXXFLAGS += -std=c++0x $(shell pkg-config --cflags jellyfish-2.0)
        LDFLAGS += $(shell pkg-config --libs jellyfish-2.0) -Wl,--rpath=$(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L//g')
    endif
endif

DEPS = kseq/kseq.h
SOURCES = jmers.cc
OBJ = $(SOURCES:.cc=.o)
OUTPUT = jmers

all: $(OUTPUT)

%.o: %.cc $(DEPS)
	$(CC) -c -o $@ $< $(CXXFLAGS)

$(OUTPUT): $(OBJ)
	$(CC) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJ)
