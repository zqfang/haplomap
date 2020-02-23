# CFLAGS=-ggdb  -m64
# CFLAGS=-ggdb -pg -O2

# CFLAGS= -ggdb -O0 -m64

CFLAGS=-Wall -ggdb -O2 -m64

# CFLAGS=-ggdb -m64

# *** Unfortunately, optimizing phmap with more than -0 seems to trigger
# a weird bug.

clean:
	rm *.o bin/eblocks bin/ghmap

ColumnReader.o:    include/ColumnReader.h  src/ColumnReader.cpp
	g++ $(CFLAGS) src/ColumnReader.cpp -c

haploLib.o:     src/haploLib.cpp include/ColumnReader.h include/dynum.h include/util.h
	g++ $(CFLAGS) src/haploLib.cpp -c

# current
ehaploblocks.o:     ehaploblocks.cpp include/ColumnReader.h include/haploLib.h include/dynum.h include/util.h
	g++ $(CFLAGS) ehaploblocks.cpp  -I./include -c

eblocks:    ehaploblocks.o ColumnReader.o haploLib.o
	g++ $(CFLAGS) ehaploblocks.o ColumnReader.o haploLib.o -l stdc++ -o bin/eblocks

# Can't optimize this more than -O because of G++ optimization bug.
# Gets different bugs with -O!
# I think I have fixed this (uninitialized var in dotVectors).
# Oops -- optimization makes new gene expression table disappear.  Valgrind doesn't see memory errors.
quantTraitMap.o:	quantTraitMap.cpp include/haploLib.h include/ColumnReader.h include/util.h
	g++ $(CFLAGS) quantTraitMap.cpp -c -L/home/fangzq/program/gsl/lib -l gsl -l gslcblas -I/home/fangzq/program/gsl/include -I./include
#	g++ -ggdb -m64 quantTraitMap.cpp -c 
#	g++ -Wall -O6 -m64 quantTraitMap.cpp -c 


ghmap:	quantTraitMap.o ColumnReader.o haploLib.o
	g++ $(CFLAGS) quantTraitMap.o ColumnReader.o haploLib.o -l stdc++ -L/home/fangzq/program/gsl/lib -l gsl -l gslcblas -o bin/ghmap

# geneNames.o:	geneNames.cpp include/ColumnReader.h include/util.h
# 	g++ $(CFLAGS) geneNames.cpp -c

# genetest:	geneNames.o ColumnReader.o
# 	g++ $(CFLAGS) geneNames.o ColumnReader.o -l stdc++  -o genetest
