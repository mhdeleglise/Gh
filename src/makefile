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
COMMANDS= testGcomb testGdelta adjust_sx_delta factorG invli2 sublong Gcompute

OBJ = erreur.o Gcomb.o Gdelta.o limpfr.o prev_prime.o GdeltaLog.o Nk.o

CXX =  g++ -I $(INC) $(LFLAGS) $(CFLAGS) $(LIBDIR)

%.o: %.cc
	$(CXX) -c   $<

all: $(OBJ) $(COMMANDS)
	echo $(UNAME)
	echo $(CFLAGS)
	echo $(LPATH)

erreur.o:erreur.cc
	$(CXX) -O2 -c erreur.cc

Gcomb.o: Gcomb.cc
	$(CXX) -O2 -c  -Wall Gcomb.cc

Gdelta.o: Gdelta.cc
	$(CXX) -O2 -c  -Wall Gdelta.cc

Nk.o: Nk.cc
	$(CXX) -O2 -c Nk.cc

GdeltaLog.o: GdeltaLog.cc
	$(CXX) -O2 -c  -Wall GdeltaLog.cc

prev_prime.o:	prev_prime.cc
	$(CXX) -O2 -c  -Wall prev_prime.cc

testGcomb:	$(OBJ) testGcomb.cc	
	$(CXX) $(CFLAGS) -o  testGcomb  testGcomb.cc $(LPATH) $(OBJ) $(LIBS)

testGdelta:	$(OBJ) testGdelta.cc
	$(CXX) $(CFLAGS) -o  testGdelta  testGdelta.cc $(LPATH) $(OBJ) $(LIBS)

Gcompute:	$(OBJ) Gcompute.cc
	$(CXX) $(CFLAGS) -o  Gcompute  Gcompute.cc $(LPATH) $(OBJ) $(LIBS) 

factorG:	$(OBJ) factorG.cc
	$(CXX) $(CFLAGS) -o  factorG  factorG.cc $(LPATH) $(OBJ) $(LIBS)

adjust_sx_delta: $(OBJ) adjust_sx_delta.cc
	$(CXX) $(CFLAGS) -o adjust_sx_delta  adjust_sx_delta.cc $(LPATH) $(OBJ) $(LIBS)

invli2:	invli2.cc $(OBJ)
	$(CXX) -o invli2 invli2.cc $(LPATH) $(OBJ) $(LIBS)

sublong:	sublong.cc $(OBJ)
	$(CXX) -o sublong sublong.cc $(LPATH) $(OBJ) $(LIBS)

limpfr.o:
	$(CXX) -c limpfr.cc

clean:
	rm -f *.o *.~  $(COMMANDS)


