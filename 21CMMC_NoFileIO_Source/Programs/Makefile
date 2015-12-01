# C compiler and flags
CPPFLAGS = -I/usr/local/include
LDFLAGS = -lgsl -lgslcblas -lfftw3f_omp -lfftw3f -lm
CC = gcc -fopenmp -Ofast
#CC = gcc -fopenmp

UNAME := $(shell uname)

# Directory hosting the non-project specific code
COSMO_DIR = ../Cosmo_c_files

# Directory containing parameter files, definint relevant constants
PARAMETER_DIR = ../Parameter_files

# non-project dependencies
COSMO_FILES = 	${PARAMETER_DIR}/COSMOLOGY.H \
	${COSMO_DIR}/ps.c \
	${COSMO_DIR}/cosmo_progs.c \
	${COSMO_DIR}/misc.c \
	${PARAMETER_DIR}/INIT_PARAMS.H \
	${PARAMETER_DIR}/ANAL_PARAMS.H \
	${PARAMETER_DIR}/HEAT_PARAMS.H \

# object files
OBJ_FILES = drive_21cmMC_streamlined \
#########################################################################

ifeq ($(UNAME), Linux)
drive_21cmMC_streamlined: drive_21cmMC_streamlined.c \
	${COSMO_FILES} \
	
	${CC} ${CPPFLAGS} -fPIC -c drive_21cmMC_streamlined.c ${LDFLAGS}
	
	${CC} -shared -Wl,-soname,libdrive_21cmMC_streamlined.so.1 -o libdrive_21cmMC_streamlined.so drive_21cmMC_streamlined.o ${LDFLAGS}
endif
ifeq ($(UNAME), Darwin)
drive_21cmMC_streamlined: drive_21cmMC_streamlined.c \
        ${COSMO_FILES} \

	${CC} ${CPPFLAGS} -o libdrive_21cmMC_streamlined.so -shared drive_21cmMC_streamlined.c ${LDFLAGS}
endif

clean:
	rm *.o *.so
