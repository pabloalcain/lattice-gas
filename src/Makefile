# -*- Makefile -*-
SHELL=/bin/sh
CC=gcc
CFLAGS=-Wall -pg -O3 -ffast-math -fPIC #-fomit-frame-pointer
LDLIBS=-lm

# list of source files
SRC=evolution.c
VPATH=src/
OBJ=evolution.o

default: liblatgas.so

clean:
	rm -f src/evolution.o

# linker rule

liblatgas.so: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ -shared $(LDLIBS)

latgas.e: $(OBJ) sim.c
	$(CC) -o $@ $(CFLAGS) $^ $(LDLIBS)
# compilation

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.c: %.h
