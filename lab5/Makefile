CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso= ..



esegui:
	make main.exe
	./main.exe
	make clean

main.exe : main.o $(percorso)/random.o lib.o $(percorso)/function.o $(percorso)/block_stat.o
	$(CC) $(percorso)/random.o main.o $(percorso)/function.o $(percorso)/block_stat.o lib.o -o main.exe

main.o : main.cpp 
	$(CC) -c main.cpp -o main.o $(CFLAGS)

$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h $(percorso)/function.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)

lib.o : lib.cpp lib.h $(percorso)/random.h
	$(CC) -c lib.cpp -o lib.o $(CFLAGS)

$(percorso)/function.o : $(percorso)/function.cpp $(percorso)/function.h  
	$(CC) -c $(percorso)/function.cpp -o $(percorso)/function.o $(CFLAGS)

$(percorso)/block_stat.o : $(percorso)/block_stat.cpp $(percorso)/block_stat.h 
	$(CC) -c $(percorso)/block_stat.cpp -o $(percorso)/block_stat.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
