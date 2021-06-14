CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: es08_1.exe es08_2.exe

es08_1.exe : es08_1.o random.o blocking.o MCMC.o
	$(CC) random.o blocking.o MCMC.o es08_1.o -o es08_1.exe
es08_1.o : es08_1.cpp
	$(CC) -c es08_1.cpp -o es08_1.o $(CFLAGS)
random.o : random.cpp random.h
es08_2.exe : es08_2.o random.o blocking.o MCMC.o
	$(CC) random.o blocking.o MCMC.o es08_2.o -o es08_2.exe
es08_2.o : es08_2.cpp
	$(CC) -c es08_2.cpp -o es08_2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
MCMC.o : MCMC.cpp random.cpp
	$(CC) -c MCMC.cpp -o MCMC.o $(CFLAGS)
blocking.o : blocking.cpp blocking.h
	$(CC) -c blocking.cpp -o blocking.o $(CFLAGS)

clean:
	rm *.o *.exe *.out
