#
# Basic Makefile for jmers. Currently it doesn't work on OS X, but that will be 
# a later implementation.
#
# This file is heavily based on Jellyfish/examples/query_per_sequence/Makefile.
#
# By: Martin Norling, 2016
#

CC = g++
CXXFLAGS = $(shell pkg-config --cflags jellyfish-2.0) -std=c++0x -Wall -O3
LDFLAGS = -lz $(shell pkg-config --libs jellyfish-2.0) -Wl,--rpath=$(shell pkg-config --libs-only-L jellyfish-2.0 | sed -e 's/-L//g')

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
