CC = g++
CFLAGS = -Wall -O3 --std=c++11

es05_1.exe : es05_1.o random.o blocking.o MCMC.o
	$(CC) random.o blocking.o MCMC.o es05_1.o -o es05_1.exe
es05_1.o : es05_1.cpp
	$(CC) -c es05_1.cpp -o es05_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
MCMC.o : MCMC.cpp MCMC.h
	$(CC) -c MCMC.cpp -o MCMC.o $(CFLAGS)
blocking.o : blocking.cpp blocking.h
	$(CC) -c blocking.cpp -o blocking.o $(CFLAGS)

clean :
	rm *.o *.exe UnifPsi100.txt GaussPsi100.txt UnifPsi210.txt GaussPsi210.txt *.out
