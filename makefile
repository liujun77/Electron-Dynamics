CC = g++
FLAG = -O3 -fopenmp -Wall -DDEBUG_TIME -DDEBUG -DOPENMP
OBJ = main.o const.o states.o readin.o save.o bessel_utilities.o\
	  legendre_rule_fast.o convolution.o
LIB = -lgsl -lgslcblas

te : $(OBJ)
	$(CC) $(FLAG) -o te $(OBJ) $(LIB)
main.o : main.cpp
	$(CC) $(FLAG) -c main.cpp
states.o : states.cpp
	$(CC) $(FLAG) -c states.cpp
convolution.o : convolution.cpp
	$(CC) $(FLAG) -c convolution.cpp
bessel_utilities.o : bessel_utilities.cpp
	$(CC) $(FLAG) -c bessel_utilities.cpp
const.o : const.cpp
	$(CC) $(FLAG) -c const.cpp
readin.o : readin.cpp
	$(CC) $(FLAG) -c readin.cpp
save.o : save.cpp
	$(CC) $(FLAG) -c save.cpp
legendre_rule_fast.o : legendre_rule_fast.cpp
	$(CC) $(FLAG) -c legendre_rule_fast.cpp

.PHONY : clean
clean:
	rm $(OBJ)
