CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso=..



esegui:
	make Gas
	make Liquid
	make Solid

all_termalized:
	make Gas_termalized
	make Liquid_termalized
	make Solid_termalized


Gas: Monte_Carlo_NVT.exe
	cp gas/input.gas input.dat
	./Monte_Carlo_NVT.exe gas
	make clean
	mv output.* gas/
	mv config.final gas/


Liquid: Monte_Carlo_NVT.exe
	cp liquid/input.liquid input.dat
	./Monte_Carlo_NVT.exe liquid
	make clean
	mv output.* liquid/
	mv config.final liquid/


Solid: Monte_Carlo_NVT.exe
	cp solid/input.solid input.dat
	./Monte_Carlo_NVT.exe solid
	make clean
	mv output.* solid/
	mv config.final solid/



Gas_termalized: Monte_Carlo_NVT.exe
	rm output.*
	cp gas/input_term.gas input.dat
	cp gas/config.final config.final
	./Monte_Carlo_NVT.exe gas
	make clean
	mv output.* gas/
	mv config.final gas/



Liquid_termalized: Monte_Carlo_NVT.exe
	rm output.*
	cp liquid/input_term.liquid input.dat
	cp liquid/config.final config.final
	./Monte_Carlo_NVT.exe liquid
	make clean
	mv output.* liquid/
	mv config.final liquid/


Solid_termalized: Monte_Carlo_NVT.exe
	rm output.*
	cp solid/input_term.solid input.dat
	cp solid/config.final config.final
	./Monte_Carlo_NVT.exe solid
	make clean
	mv output.* solid/
	mv config.final solid/




Monte_Carlo_NVT.exe : Monte_Carlo_NVT.o $(percorso)/random.o
	$(CC) $(percorso)/random.o Monte_Carlo_NVT.o -o Monte_Carlo_NVT.exe
Monte_Carlo_NVT.o : Monte_Carlo_NVT.cpp
	$(CC) -c Monte_Carlo_NVT.cpp -o Monte_Carlo_NVT.o $(CFLAGS)
$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_NVT.exe
