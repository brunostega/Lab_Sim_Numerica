CC = mpicxx
GG = g++
CFLAGS = -Wall -O3 --std=c++11
percorso= ../..



esegui: 
	make Circle
	make Square

Circle: main.exe
	cp input.circle input.dat
	mpiexec -np 6 ./main.exe
	mv *.circle.* Circle/	
	make clean

Square: main.exe
	cp input.square input.dat
	mpiexec -np 6 ./main.exe
	mv *.square.* Square/	
	make clean

main.exe: main.o $(percorso)/random.o 
	$(CC)  -o main.exe $(percorso)/random.o main.o $(percorso)/function.o $(percorso)/block_stat.o

main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h $(percorso)/function.h
	$(GG) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)


clean:
	rm *.o *.exe
