CC	= clang++
SRC	= lat_fact.cpp
CFLAGS	= -Ofast -fopenmp -mtune=native -march=native -pg -std=c++2b
LDFLAGS = -lntl

all:
	${CC} ${CFLAGS} ${SRC} ${LDFLAGS}
