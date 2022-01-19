CC = g++
CFLAGS = -Wall -O3 --std=c++11
percorso= ../..


all:
	make Gas
	make Liquid
	make Solid

all_termalized:
	make Gas_termalized
	make Liquid_termalized
	make Solid_termalized


Gas: main.exe
	rm output.*
	rm gas/output_temp.dat
	touch gas/output_temp.dat
	cp gas/input_gas.dat input.dat
	./main.exe gas
	cp output.* gas/
	make clean


Liquid: main.exe
	rm output.*
	rm liquid/output_temp.dat
	touch liquid/output_temp.dat
	cp liquid/input_liquid.dat input.dat
	./main.exe liquid
	cp output.* liquid/
	make clean


Solid: main.exe
	rm output.*
	rm solid/output_temp.dat
	touch solid/output_temp.dat
	cp solid/input_solid.dat input.dat
	./main.exe solid
	cp output.* solid/
	make clean


Gas_termalized: 
	cp gas/input_gas_term.dat input.dat
	./main.exe gas
	make clean

Liquid_termalized: 
	cp liquid/input_liquid_term.dat input.dat
	./main.exe liquid
	make clean

Solid_termalized:
	cp solid/input_solid_term.dat input.dat
	./main.exe solid
	make clean

normal: 
	rm output_temp.dat
	touch output_temp.dat
	./main.exe

main.exe : main.o $(percorso)/random.o $(percorso)/function.o $(percorso)/block_stat.o
	$(CC) $(percorso)/random.o main.o $(percorso)/function.o $(percorso)/block_stat.o -o main.exe

main.o : main.cpp main.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)

$(percorso)/random.o : $(percorso)/random.cpp $(percorso)/random.h $(percorso)/function.h
	$(CC) -c $(percorso)/random.cpp -o $(percorso)/random.o $(CFLAGS)

$(percorso)/function.o : $(percorso)/function.cpp $(percorso)/function.h  
	$(CC) -c $(percorso)/function.cpp -o $(percorso)/function.o $(CFLAGS)

$(percorso)/block_stat.o : $(percorso)/block_stat.cpp $(percorso)/block_stat.h 
	$(CC) -c $(percorso)/block_stat.cpp -o $(percorso)/block_stat.o $(CFLAGS)

clean :
	rm *.o main.exe 
