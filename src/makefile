INC= ~deleglis/Pix2015/Include
LIBDIR=-L ~deleglis/Pix2015
LIBS =  -lgmpxx -lgmp -lmpfr -lsieve
CFLAGS = -O2  -Wall

UNAME :=$(shell uname)

ifeq ($(UNAME),Darwin)
	CFLAGS+=-stdlib=libstdc++
	LPATH=-L/usr/local/lib
else
	LPATH=
endif

SOURCES= erreur.cc primes1000000.cc Gcomb.cc Gdelta.cc
COMMANDS= testGcomb testGdelta adjust_sx_delta factorG invli2 sublong\
	Gcompute left_prime

OBJ = erreur.o Gcomb.o Gdelta.o limpfr.o prev_prime.o GdeltaLog.o Nk.o\
	psi.o psi_theta_diff.o  theta.o

CXX =  g++ -I $(INC) $(LFLAGS) $(CFLAGS)

%.o: %.cc
	$(CXX) -c   $<

all: $(OBJ) $(COMMANDS)
	echo $(UNAME)
	echo $(CFLAGS)
	echo $(LPATH)

left_prime:
	$(CXX) -o left_prime left_prime.cc -lgmp

testGcomb:	$(OBJ) testGcomb.cc	
	$(CXX) $(CFLAGS) -o  testGcomb  testGcomb.cc $(LPATH) $(OBJ) $(LIBS) $(LIBDIR)

testGdelta:	$(OBJ) testGdelta.cc
	$(CXX) $(CFLAGS) -o  testGdelta  testGdelta.cc $(LPATH) $(OBJ) $(LIBS)  $(LIBDIR)

Gcompute:	$(OBJ) Gcompute.cc
	$(CXX) $(CFLAGS) -o  Gcompute  Gcompute.cc $(LPATH) $(OBJ) $(LIBS)  $(LIBDIR)

factorG:	$(OBJ) factorG.cc
	$(CXX) $(CFLAGS) -o  factorG  factorG.cc $(LPATH) $(OBJ) $(LIBS)  $(LIBDIR)

adjust_sx_delta: $(OBJ) adjust_sx_delta.cc
	$(CXX) $(CFLAGS) -o adjust_sx_delta  adjust_sx_delta.cc $(LPATH) $(OBJ) $(LIBS)  $(LIBDIR)

invli2:	invli2.cc $(OBJ)
	$(CXX) -o invli2 invli2.cc $(LPATH) $(OBJ) $(LIBS)  $(LIBDIR)

sublong:	sublong.cc $(OBJ)
	$(CXX) -o sublong sublong.cc $(LPATH) $(OBJ) $(LIBS)    $(LIBDIR)

limpfr.o:
	$(CXX) -c limpfr.cc

clean:
	rm -f *.o *.~  $(COMMANDS)


