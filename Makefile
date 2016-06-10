all:
	cd sieve2015 && $(MAKE)
	cd src && $(MAKE)


clean:
	cd src && $(MAKE) clean
	cd sieve2015 && $(MAKE) clean
