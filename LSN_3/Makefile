CC = g++
CFLAGS = -Wall -O3 --std=c++11

es03_1.exe : es03_1.o random.o brownian_motion.o blocking.o
	$(CC) random.o brownian_motion.o blocking.o es03_1.o -o es03_1.exe
es03_1.o : es03_1.cpp
	$(CC) -c es03_1.cpp -o es03_1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
brownian_motion.o : brownian_motion.cpp brownian_motion.h
	$(CC) -c brownian_motion.cpp -o brownian_motion.o $(CFLAGS)
blocking.o : blocking.cpp blocking.h
	$(CC) -c blocking.cpp -o blocking.o $(CFLAGS)
clean :
	rm *.o *.exe *.txt *.out
