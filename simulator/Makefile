GPP = g++
OPT = -O4 -g
OBJ = main.o RiboModel.o Ribosome.o CodonTranslator.o
VPATH =../lib
LIBPATH = $(VPATH)

TASEPsimulator.exe: $(OBJ)
	$(GPP) $(OPT) $(OBJ) -o $@

main.o: main.cpp
	$(GPP) $(OPT) -c $< -o $@

RiboModel.o: $(LIBPATH)/RiboModel.cpp
	$(GPP) $(OPT) -c $< -o $@

Ribosome.o: $(LIBPATH)/Ribosome.cpp 
	$(GPP) $(OPT) -c $< -o $@

CodonTranslator.o: $(LIBPATH)/CodonTranslator.cpp 
	$(GPP) $(OPT) -c $< -o $@

clean:
	rm *.o
	rm *~
