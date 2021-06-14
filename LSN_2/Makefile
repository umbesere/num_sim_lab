CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: es02_1.exe es02_2.exe

es02_2.exe : es02_2.o random.o random_walk.o
	$(CC) random.o random_walk.o es02_2.o -o es02_2.exe
es02_2.o : es02_2.cpp
	$(CC) -c es02_2.cpp -o es02_2.o $(CFLAGS)
es02_1.exe : es02_1.o random.o
	$(CC) random.o es02_1.o -o es02_1.exe
es02_1.o : es02_1.cpp
	$(CC) -c es02_1.cpp -o es02_1.o $(CFLAGS)
random_walk.o : random_walk.cpp random_walk.h
	$(CC) -c random_walk.cpp -o random_walk.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt seed.out
exe:
	./es02_1.exe ./es02_2.exe
