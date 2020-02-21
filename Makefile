# CFLAGS=-ggdb  -m64
# CFLAGS=-ggdb -pg -O2

# CFLAGS= -ggdb -O0 -m64

CFLAGS=-Wall -ggdb -O2 -m64

# CFLAGS=-ggdb -m64

# *** Unfortunately, optimizing phmap with more than -0 seems to trigger
# a weird bug.

#clean:    
#	rm *.o eblocks ghmap

ColumnReader.o:    ColumnReader.h  ColumnReader.cpp
	g++ $(CFLAGS) ColumnReader.cpp -c

haploLib.o:     haploLib.cpp ColumnReader.h dynum.h util.h
	g++ $(CFLAGS) haploLib.cpp -c

# current
ehaploblocks.o:     ehaploblocks.cpp ColumnReader.h haploLib.h dynum.h util.h
	g++ $(CFLAGS) ehaploblocks.cpp -c

eblocks:	ehaploblocks.o ColumnReader.o haploLib.o
	g++ $(CFLAGS) ehaploblocks.o ColumnReader.o haploLib.o -l stdc++ -o eblocks

# Can't optimize this more than -O because of G++ optimization bug.
# Gets different bugs with -O!
# I think I have fixed this (uninitialized var in dotVectors).
# Oops -- optimization makes new gene expression table disappear.  Valgrind doesn't see memory errors.
quantTraitMap.o:	quantTraitMap.cpp haploLib.h ColumnReader.h util.h
	g++ $(CFLAGS) quantTraitMap.cpp -c 
#	g++ -ggdb -m64 quantTraitMap.cpp -c 
#	g++ -Wall -O6 -m64 quantTraitMap.cpp -c 


ghmap:	quantTraitMap.o ColumnReader.o haploLib.o
	g++ $(CFLAGS) quantTraitMap.o ColumnReader.o haploLib.o -l stdc++ -L/usr/lib -l gsl -l gslcblas -o ghmap

geneNames.o:	geneNames.cpp ColumnReader.h util.h
	g++ $(CFLAGS) geneNames.cpp -c

genetest:	geneNames.o ColumnReader.o
	g++ $(CFLAGS) geneNames.o ColumnReader.o -l stdc++  -o genetest
