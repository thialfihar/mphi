CC=g++
CPPFLAGS=-std=c++11 -Ofast -fopenmp

all: find_lambda map

functions.o: functions.cpp functions.h
	$(CC) $(CPPFLAGS) -c -o functions.o functions.cpp

benchmark.o: benchmark.cpp Matrix.h GeneralMatrix.h Mutex.h
	$(CC) $(CPPFLAGS) -c -o benchmark.o benchmark.cpp

find_lambda.o: find_lambda.cpp Matrix.h Mutex.h
	$(CC) $(CPPFLAGS) -c -o find_lambda.o find_lambda.cpp

map.o: map.cpp Matrix.h Mutex.h
	$(CC) $(CPPFLAGS) -c -o map.o map.cpp

find_lambda: functions.o find_lambda.o
	$(CC) $(CPPFLAGS) -o find_lambda functions.o find_lambda.o -lgmp

benchmark: functions.o benchmark.o
	$(CC) $(CPPFLAGS) -o benchmark functions.o benchmark.o -lgmp

map: functions.o map.o
	$(CC) $(CPPFLAGS) -o map functions.o map.o -lgmp

clean:
	rm *.o
