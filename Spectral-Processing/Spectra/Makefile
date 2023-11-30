CC=gcc
CFLAGS=-I.
# Uncomment the following line to enable OpenMP
# OMPFLAG=-fopenmp
DEPS=atmos.h constant.h include.h input.h nrutil.h opac.h two_stream.h
OBJ=geometry.o interpol.o main_rt_no_scat.o planck.o read_t_p_doppler.o readchemtable.o readopactable.o totalopac.o rt_emission_2stream_multiphase.o utils.o nrutil.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OMPFLAG)

rt_emission_aerosols.exe: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(OMPFLAG) -lm

clean:
	rm -f rt_emission_aerosols.exe *.o
