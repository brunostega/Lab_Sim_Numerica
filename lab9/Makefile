CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso=..


all: 
	make circle 
	make square

circle: main.exe
	./main.exe circle
	make clean

square: main.exe
	./main.exe square
	make clean

main.exe : main.o $(percorso)/random.o  $(percorso)/function.o $(percorso)/block_stat.o
	$(CC) $(percorso)/random.o main.o $(percorso)/function.o $(percorso)/block_stat.o  -o main.exe

main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h $(percorso)/function.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)


$(percorso)/function.o : $(percorso)/function.cpp $(percorso)/function.h  
	$(CC) -c $(percorso)/function.cpp -o $(percorso)/function.o $(CFLAGS)

$(percorso)/block_stat.o : $(percorso)/block_stat.cpp $(percorso)/block_stat.h 
	$(CC) -c $(percorso)/block_stat.cpp -o $(percorso)/block_stat.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
