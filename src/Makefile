# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -pg -O3 -ffast-math -fPIC -fopenmp
LDLIBS=-lm -lgomp

# list of source files
SRC=evolution.c
VPATH=src/
OBJ=evolution.o

default: liblatgas.so

clean:
	rm -f evolution.o latgas.e liblatgas.so

# linker rule

liblatgas.so: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ -shared $(LDLIBS)

latgas.e: $(OBJ) sim.c
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)
# compilation

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.c: %.h
