UNAME := $(shell uname)
ifeq ($(UNAME),Darwin)
	LFLAGS+=-stdlib=libstdc++
	LPATH=-L/opt/local/lib
endif

include dflags.mk

INC=-I../Include
LPATH=-L../  -lgmp

CXX=g++ $(DFLAGS) -O3 -Wall $(INC) 

OBJS=    arith_long128.o arith_functions.o\
	bit_table.o bit_table_cnte.o \
	sieve_access_frame.o presieved_primes.o sieve128_pg.o\
	sieve_large.o nth_rightp.o nth_leftp.o\
	primes.o primes6.o nextPrime.o

%	: %.cc .headers
	$(CXX) $< $(OBJS) $(LFLAGS) $(LPATH) -o $@

TEST_SRC:=$(shell ls tst*.cc)
TEST=$(patsubst %.cc, %, $(TEST_SRC))


libsieve.a:     $(OBJS)   
	ar -r libsieve.a $(OBJS)
	mv libsieve.a ../
	make tests

tests: ../libsieve.a $(TEST)


# headers files
.headers: ../Include/arith_long128.h ../Include/arith_64_32.h ../Include/bit_table.h ../Include/bit_table_cnte.h\
	../Include/sieve_access_frame.h ../Include/arith_basic_functions.h\
	../Include/arith_structured_types.h ../Include/sieve_slices.h \
	../Include/mylib.h ../Include/interval_sieve.h ../Include/interval_sieve.hpp\
	../Include/prime_generator.h ../Include/sieve_slices.hpp
	@touch $@


clean:
	rm -f *.o $(TEST)

# have a look at https://www.gnu.org/software/make/manual/html_node/index.html#SEC_Contents


