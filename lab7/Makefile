CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso=..

esegui: Monte_Carlo_NVT.exe
	rm output.*
	./Monte_Carlo_NVT.exe
	make clean

all:
	make Gas
	make Liquid
	make Solid

all_termalized:
	make Gas_termalized
	make Liquid_termalized
	make Solid_termalized


Gas: Monte_Carlo_NVT.exe
	rm output.*
	cp gas/input.gas input.dat
	./Monte_Carlo_NVT.exe gas
	make clean
	cp output.* gas/
	cp config.final gas/



Liquid: Monte_Carlo_NVT.exe
	rm output.*
	cp liquid/input.liquid input.dat
	./Monte_Carlo_NVT.exe liquid
	make clean
	cp output.* liquid/
	cp config.final liquid/


Solid: Monte_Carlo_NVT.exe
	rm output.*
	cp solid/input.solid input.dat
	./Monte_Carlo_NVT.exe solid
	make clean
	cp output.* solid/
	cp config.final solid/



Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o $(percorso)/random.o
	$(CC) $(percorso)/random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe
