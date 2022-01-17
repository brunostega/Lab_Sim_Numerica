CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso= ..

Metropolis: Monte_Carlo_ISING_1D.exe
	
	rm outputMetropolis/output.*
	./Monte_Carlo_ISING_1D.exe metro
	make clean
	cp ../seed.out outputMetropolis/seed.out

Gibbs: Monte_Carlo_ISING_1D.exe

	rm outputGibbs/output.*
	./Monte_Carlo_ISING_1D.exe gibbs
	make clean
	cp ../seed.out outputGibbs/seed.out

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o $(percorso)/random.o
	$(CC) $(percorso)/random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe 
