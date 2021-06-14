CC = g++
CFLAGS = -Wall -O3 --std=c++11

es9.exe : es9.o random.o genetic.o
	$(CC) random.o genetic.o es9.o -o es9.exe $(CFLAGS)
es9.o : es9.cpp
	$(CC) -c es9.cpp -o es9.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
genetic.o : genetic.cpp genetic.h
	$(CC) -c genetic.cpp -o genetic.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt *.out
