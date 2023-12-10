# Spectral-Processing

*************************************************
************       USER GUIDE        ************
*************************************************
This is a guide for how to use the files here

Recently, Isaac and Alex have changed how to run the post processing. It is now much easier and faster. 


*************************************************
*********         Before Running      **********
*************************************************

The code needs some data files to run that are too big to keep in github. Isaac put these on Turbo
Running download_all_data_files_from_turbo.py should get you the files
Look at the bottom of the file in files_to_download to see which ones its grabbing

Once the files are downloaded, move the SET_X of data to Spectra/DATA, and move the SCATTERING_DATA to Spectra/SCATTERING_DATA

Also, make sure to put all the GCMs that you want to run in GCM-OUTPUT

*************************************************
************         HOW TO RUN      ************
*************************************************

Spectral-Processing/run_entire_suite.py
    This will run the entire suite

    The entire post processing suite will run by typing 'sbatch Run_all_sbatch' in your command line. A good go to for runtime is at least 24 hours. 
    The benefit is that it is very fast to get everything in parallel, and now everything is run by one click of the button. 
    
    Be sure to specify which phases you want run in run_entire_suite.py . 
	
    There are three parts to the post processing —altitude regridding, init files, and the final calculations. 
    STEP_ONE: Altitude regridding takes around an hour and submits one job. It creates 4 files for every planet and places them in /PLANET_MODELS/ . 
    STEP_TWO: The init files are placed in /Spectra/DATA/ . It will create a separate init file for every phase you run. These take around 20 minutes. 
    STEP_THREE: Finally, the ending calculations will be placed in /Spectra/OUT/ . There will be two data files for each phase ran. These can take 12 hours depending on the cloudiness of your model. 
    All three above steps are run at the same time. The code will automatically determine if the altitude regridding and the init files need to be made. This means that if you want to redo them, be sure to delete / move the existing files out. It will not overwrite files in /PLANET_MODELS/ or /Spectra/DATA/ . 

    The main part of the program that this is calling is Spectra/run_spectra.py
    This will run all the subprograms for regridding, interpolating, and finally running the post-processing

    Some things to keep in mind
        NLAT: The final number of layers the model will have, try 250 or 500
        opacity_files: Which set of opacity files you're using

    The code will grab all necessary information from the fort.7 . 

    Also, if you care about getting the cloud dependent wavelength properties, you do have to set that in run_spectra.py
    This is only for specific visualization stuff, not the general postprocessing


Visualizations/create_all_figures.py
    This will create all the figures that I've so far coded up

After all of this the empty folders will probably need to be deleted.

*************************************************
************         FOLDERS        *************
*************************************************

First, there are several folders:
- Figures            || The figures
- Visualizations     || The code for making the figures
- FINISHED_SPECTRA   || This is where all the finished output spectra should go
- GCM-OUTPUT         || This is where all the GCM outputs from the RM-GCM should go
- PLANET_MODELS      || This is where all the regridded GCMS and interpolations go
- Spectra            || This is the main folder where the spectral procressing happens. DON'T PUT DATA HERE


*************************************************
************     FILE NAMING       **************
*************************************************

The new version of the post processing grabs all necessary information from the fort.7 . This means that specific naming conventions for the GCM-OUTPUT subfolders are no longer necessary. 

However, it is important that planets are titled differently to prevent overwriting data output files. 

Most importantly, check that the data output that is printed to the terminal matches what you want!
It will print out all the important characterists!


*************************************************
********      OPACITY VERSIONS              *****
*************************************************

This code can work for a number of different resolution and temperature regimes.
However, this requires swapping out the EOS and opacity files, as well as changing the code itself to read in these.
I don't think that this is a good idea, but it would be a huge amount of work to change it.
In the meantime, the code will simply swap out the files that need to be changed.
These files are in OPAC_CODE_VERSIONS and the different opac and EOS files are in DATA

The files are broken up into the following sets:
SET_1: The low temp files for GJ1214b 
SET_2: The files for HD189 + HD209. These are pretty standard hot jupiter files
SET_3: The high res files, identical to SET_2 but high res

AS YOU ADD STUFF TO THESE SETS PLEASE UPDATE THIS README!!!!

*************************************************
************          NOTES         *************
*************************************************

In /Spectra/ , there will be (6 * number of phases + 2) number of files created. This is normal and is required for running the phases in parallel. You can easily delete them by running Clean_suite.py (by typing 'python3 Clean_suite.py' in command line). There will also be many slurm files created. Sorry. 

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
