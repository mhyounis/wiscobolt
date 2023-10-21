# Compiler settings
FC = gfortran
FCFLAGS = -fopenmp -Wl,--stack,100000000 -fdefault-real-8

# Modules
MODULES = math.o physics.o 1_user_input.o 2_mesh.o 3_shape_functions.o \
          4_angular_discretization.o 5_sweep_order.o 6_spatial_inner_products.o \
          7_energy_discretization.o MGXS.o FEXS.o 8_sources.o 8_ele_sources.o 9_iteration.o \
          9_t1_ele_iteration.o 9_t2_ele_iteration.o post_processing.o

# Program
PROGRAM = wiscobolt

all: $(PROGRAM)

$(PROGRAM): $(MODULES) solver.o
	$(FC) $(FCFLAGS) $(MODULES) solver.o -o $(PROGRAM)

%.o: %.f08
	@$(FC) $(FCFLAGS) -c $< -o $@ > NUL

run: $(PROGRAM)
	@./$(PROGRAM) $(input)

clean:
	@del /Q $(PROGRAM) $(PROGRAM).exe *.o *.mod