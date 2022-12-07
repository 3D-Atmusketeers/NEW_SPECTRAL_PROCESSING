# Spectral-Processing
This is a guide for how to use the files here


Recently, Isaac has added a lot of stuff to the files. Hopefully this makes it way easier to use, but it also may cause problems if not documented.

First, there are several folders:
- FINISHED_SPECTRA   || This is where all the finished output spectra should go
- GCM-OUTPUT         || This is where all the GCM outputs from the RM-GCM should go
- PLANET_MODELS      || This is where all the regridded GCMS and interpolations go
- Spectra            || This is the main folder where the spectral procressing happens. DON'T PUT DATA FILES HERE


One thing is that these files now expect a very specific format for the GCM names. This allows for super efficiency, but you have to be careful.

NORMAL PLANET NAME:
PLANETNAME_CLOUDTYPE_CLOUDNUMBER_HAZETYPE

So you would name it something like GJ1214b_NUC_25LAYERS_SOOT_DENSE

The code will automatically parse from this that the planet has nucleation limited clouds, 25 cloud layers, and soot hazes
If you don't want any hazes make sure to put clear (or edit the run_spectra file)

Use "ALL" for all clouds, "NUC" for nucleation, and don't include the word "CLOUDS" if you don't want clouds
Use "SOOT" for soot hazes, "THOLIN" for tholin hazes, and "CLEAR" for no hazes.
Use the keyword "DENSE" to change the MTLX to 1.0. It will assume 0.1 otherwise.

Also be very careful if you put different metallities that the model can confuse metallicity and cloud layers.

Most importantly, check that the data output that is printed to the terminal matches what you want!
It will print out all the important characterists!

The code currently can NOT do more than one planet system at a time. This is because star radius and effective temp and
stuff like that needs to be changed per run.

The run_spectra will allow you to run planets on a single core.
However, if you want to run lots of planets at once, use set_up_spectra_folders.py
This will run everything at once. It is very powerful, but will copy over huge amounts of data, so be careful.

Files that are needed:
fort.26, fort.50, fort.7, fort.62


IMPORTANT NOTE:
set_up_spectra_folders WILL copy over folders and stuff, but you will have to reorganize once everything is run
that way everything goes to the right place