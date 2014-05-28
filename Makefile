CC=g++
CPPFLAGS=-std=c++11 -Ofast

all: benchmark find_lambda

functions.o: functions.cpp functions.h
	$(CC) $(CPPFLAGS) -c -o functions.o functions.cpp

benchmark.o: benchmark.cpp Matrix.h GeneralMatrix.h
	$(CC) $(CPPFLAGS) -c -o benchmark.o benchmark.cpp

find_lambda.o: find_lambda.cpp Matrix.h
	$(CC) $(CPPFLAGS) -c -o find_lambda.o find_lambda.cpp

find_lambda: find_lambda.o benchmark.o
	$(CC) $(CPPFLAGS) -o find_lambda functions.o find_lambda.o -lgmp

benchmark: functions.o benchmark.o
	$(CC) $(CPPFLAGS) -o benchmark functions.o benchmark.o -lgmp

clean:
	rm *.o
