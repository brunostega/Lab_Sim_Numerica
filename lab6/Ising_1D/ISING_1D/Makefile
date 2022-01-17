CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso= /home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi

esegui: Monte_Carlo_ISING_1D.exe
	rm E_final.dat
	touch E_final.dat
	rm Heat_final.dat
	touch Heat_final.dat
	rm M_final.dat
	touch M_final.dat
	rm X_final.dat
	touch X_final.dat
	./Monte_Carlo_ISING_1D.exe

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o $(percorso)/random.o
	$(CC) $(percorso)/random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe
