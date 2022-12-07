## Caleb Harada's version of the code.

This version adds aerosol scattering to high-res spectra calculation. 

Edit the input header file:

-name output files

-specify the T-P profile with aerosols properties

-specify orbital phases

-switch Doppler shifts on or off

-switch clouds on or off

### Compiling:

`$ make rt_emission_aerosols.exe`

-OR-

`$ cc geometry.c interpol.c main_rt_no_scat.c planck.c readchemtable.c readopactable.c totalopac.c rt_emission_unshift_lon_fixphase.c utils.c nrutil.c read_t_p_doppler.c -lm -ansi`

### To run the code:

`$ ./rt_emission_aerosols.exe`

-OR-

`$ ./a.out`

### Clean working directory:

`$ make clean`
