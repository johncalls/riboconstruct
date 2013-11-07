################################################################################
#
#  Modify this part and set the path of your Vienna directory
#  The Vienna RNA Package can be downloaded at:
#                       http://www.tbi.univie.ac.at/~ivo/RNA/
#  and has to be configured (and installed).

VIENNA_RNA      = /usr/local
RNA_F           = /home/john/Desktop/riboconstruct/rnaf
INSTALL_DIR     = /home/john/Desktop/riboconstruct/riboconstruct/rna_f

################################################################################

CXX             = g++
CC              = gcc

CCFLAGS         = -O4 -Wall -fPIC -fopenmp
LIBFLAGS        = -shared -O4 -Wall -fopenmp

INCLUDES        = -I$(VIENNA_RNA)/include/ViennaRNA -I$(RNA_F)
LIBS            = -L$(RNA_F) -lRNAf -L$(VIENNA_RNA)/lib -lRNA
RUNPATH         = -Wl,-rpath=$(RNA_F)

################################################################################

NAME = RNAf_wrap

%.o: %.c
	$(CXX) -c $(CCFLAGS) $(INCLUDES) -o $@ $<

################################################################################

all: lib$(NAME).so

lib$(NAME).so: wrap.o
	$(CC) $(LIBFLAGS) -o $@ wrap.o $(RUNPATH) $(LIBS)

clean:
	$(RM) wrap.o lib$(NAME).so *~ gmon.out .gdb_history core
