CC	= g++

libfact_c:
	${CC} -shared -fPIC -O3 -mavx2 -fopenmp -mtune=native -march=native -mfpmath=both -unroll-loops -o libfact.so src/lat_fact.cpp -lntl
