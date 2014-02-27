# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -g -fbounds-check -Wall
#FCFLAGS = -O3

PROJDIR = .
MDIR = ${PROJDIR}/modules
OBJ = ${PROJDIR}/obj

# List of executables to be built within the package
PROGRAMS = ${MDIR}/eispack.o ${MDIR}/cornelius.o ${MDIR}/gauss_legendre.o ${MDIR}/Histogram_module.o ${MDIR}/read_f14f15.o ${MDIR}/T_from_EoS.o ${MDIR}/CF_int.o ${MDIR}/Land_Eck.o Get_part_lines Analyze_part_lines

all: $(PROGRAMS)

Get_part_lines: Get_part_lines.o read_f14f15.o
	$(FC) $(FCFLAGS) -o $@ $^ -I${OBJ}

Analyze_part_lines: Analyze_part_lines.o read_f14f15.o cornelius.o eispack.o Land_Eck.o
	$(FC) $(FCFLAGS) -o $@ $^ -I${OBJ}

%: %.o
	$(FC) $(FCFLAGS) -o $@ $^ -I${OBJ}

%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -J${OBJ} -I${OBJ}

clean:
	rm -f *.o *.mod *.MOD



