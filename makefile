compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -fsanitize=address -g -pg

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: faultCost.c utilities/readGraph6.c utilities/bitset.h utilities/hamiltonicityMethods.h
	$(compiler) -DUSE_64_BIT -o faultCost faultCost.c utilities/readGraph6.c utilities/hamiltonicityMethods.c $(flags)

# There are two different implementations of the 128-bit version. The array version generally performs faster.
128bit: faultCost.c utilities/readGraph6.c utilities/bitset.h utilities/hamiltonicityMethods.h
	$(compiler) -DUSE_128_BIT -o faultCost-128 faultCost.c utilities/readGraph6.c utilities/hamiltonicityMethods.c $(flags)

128bitarray: faultCost.c utilities/readGraph6.c utilities/bitset.h utilities/hamiltonicityMethods.h
	$(compiler) -DUSE_128_BIT_ARRAY -o faultCost-128a faultCost.c utilities/readGraph6.c utilities/hamiltonicityMethods.c $(flags)	

profile: faultCost.c utilities/readGraph6.c utilities/bitset.h utilities/hamiltonicityMethods.h
	$(compiler) -DUSE_64_BIT -o faultCost-pr faultCost.c utilities/readGraph6.c utilities/hamiltonicityMethods.c $(profileflags)

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f faultCost faultCost-128 faultCost-128a faultCost-pr

