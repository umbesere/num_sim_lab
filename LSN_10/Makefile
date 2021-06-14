CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: es10_1.exe es10_2.exe

es10_1.exe : es10_1.o random.o genetic.o
	$(CC) random.o genetic.o es10_1.o -o es10_1.exe $(CFLAGS)
es10_1.o : es10_1.cpp
	$(CC) -c es10_1.cpp -o es10_1.o $(CFLAGS)
es10_2.exe : es10_2.cpp random.h random.cpp genetic.h genetic.cpp
	mpicxx random.h random.cpp genetic.h genetic.cpp es10_2.cpp -o es10_2.exe $(CFLAGS)
try.exe : try.cpp random.h random.cpp genetic.h genetic.cpp
	mpicxx random.h random.cpp genetic.h genetic.cpp try.cpp -o try.exe $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt *.out
