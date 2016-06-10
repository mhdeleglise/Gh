all:
	cd src && $(MAKE)
	cd sieve2015 && $(MAKE)

clean:
	cd src && $(MAKE) clean
