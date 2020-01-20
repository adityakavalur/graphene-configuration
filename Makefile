CC = g++
CFLAGS = -Wall
DEPS = graphene_initialization.h graphene_nucleation.h graphene_grain_growth.h graphene_atom_generation.h graphene_grain_area_calc.h graphene_add_layer.h
OBJ = graphene_config_generator.o graphene_initialization.o graphene_nucleation.o graphene_grain_growth.o graphene_atom_generation.o graphene_grain_area_calc.o graphene_add_layer.o
EXECUTABLE = poly

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXECUTABLE): $(OBJ)
	$(CC) $(OBJ) $(CFLAGS) -o $@

clean:
	/bin/rm *.o
clear:  
	/bin/rm *.o

