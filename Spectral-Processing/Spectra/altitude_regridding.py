from tqdm import tqdm
import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter
import pandas as pd

### ----- INPUTS AND OUTPUTS ----- ###


def regrid_gcm_to_constant_alt(path, CLOUDS, planet_name, NLAT, NLON, NTAU, NLON_new, NTAU_new, HAZES, max_pressure_bar, smoothing):
    old_file = '../PLANET_MODELS/' + planet_name + '_with_clouds.txt'
    new_file = '../PLANET_MODELS/' + planet_name + '_with_clouds_regridded.txt'
    NPARAMS = 51

    #########################################
    #            END USER INPUTS                #
    #########################################

    ### ----- FIND NEW ALTITUDE GRID ----- ###
    def altitudes(data):

        # data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]


        # intial min and max altitudes

        z_min = 0

        z_max = 0


        # find highest and lowest altitude
        print ("Find the highest and lowest altitudes")

        for i in range(NLAT):

            for j in range(NLON):

                if data[i][j][0][3] > z_max:

                    z_max = data[i][j][0][3]

                if data[i][j][NTAU - 1][3] > z_min:

                    z_min = data[i][j][NTAU - 1][3]


        # define new grid of altitudes

        z_grid = np.flip(np.linspace(z_min, z_max, NTAU_new), axis=0)


        # return grid of new altitudes

        return z_grid







    ### ----- FUNCTION TO DOUBLE ALL DATA ----- ###


    def double_lons(data):
        data2 = np.copy(data)

        data3 = []


        # add 360 to copy of longitudes (does not duplicate 360)

        for row in data2:

            row[1] += 360.

            if row[1] == 360.:

                last_row = np.copy(row)

                last_row[1] += 360.

                data3.append(last_row)


        # combine copies of data        

        double_data = np.concatenate((data, data2, data3))


        # sort double_data by latitude

        double_data = double_data[double_data[:,0].argsort()]


        # sort data to match original format (sorry, this is messy and probably not efficient)
        print("Regridding to new number of layers")
        for i in tqdm(range(NLAT)):

            chunk = double_data[int(i * len(double_data) / NLAT) : int((i + 1) * len(double_data) / NLAT)]

            chunk = chunk[chunk[:,1].argsort()]

            for j in range(NLON_new):

                sub_chunk = chunk[int(j * len(chunk) / NLON_new) : int((j +1 ) * len(chunk) / NLON_new)]

                sub_chunk = sub_chunk[sub_chunk[:,2].argsort()]

                chunk[int(j * len(chunk) / (NLON_new)) : int((j + 1) * len(chunk) / (NLON_new))] = sub_chunk

            double_data[int(i * len(double_data) / NLAT) : int((i + 1) * len(double_data) / NLAT)] = chunk


        # convert bars to pascals

        #for i, row in enumerate(double_data):
        #    double_data[i][4] *= 1e+5


        # return doubled data grid

        return double_data


    ### ----- LINEAR INTERPOLATION OVER ENTIRE GRID ----- ###


    def LInterp_1d(data, data_new, z_new, param_col):

        # data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]
        # z_grid is array of length NTAU containing new altitude grid points
        """
        if param_col == 5:
            print ('here')
            print (data[0][0][:,param_col])

            param_interp = interpolate.interp1d(data[0][0][:,3], data[0][0][:,param_col], kind="linear", bounds_error=False, fill_value=0)

            param_new = param_interp(z_new)

            print (param_new)
        """


        for i in range(NLAT):

            for j in range(NLON):


                # old altitude grid at this lat, lon

                z_old = data[i][j][:,3]

                # parameter values on old altitude grid

                param_old = data[i][j][:,param_col]



                # linear interpolation function created from old altitudes and parameter values

                param_interp = interpolate.interp1d(z_old, param_old, kind="linear", bounds_error=False, fill_value=0)



                # apply interpolation function (values not on new grid set to zero)

                param_new = param_interp(z_new)



                # change parameter values in data array to the new interpolated values
                for k in range(NTAU_new):
                    if (param_col == 9 or param_col == 12 or param_col == 15 or param_col == 18):
                        data_new[i][j][k][param_col] = param_new[k] #/ (NTAU_new / NTAU)# optical depth
                    else:
                        data_new[i][j][k][param_col] = param_new[k]

        return data_new








    ### ----- LOGARITHMIC INTERPOLATION OVER ENTIRE GRID ----- ###


    def LogInterp_1d(data, data_new, z_new, param_col, integrate=False):

        # data must be an array of dimensions [NLAT x NLON x NTAU x NPARAMS]
        # z_grid is array of length NTAU containing new altitude grid points




        if integrate == True:


            # integrate

            for i in range(NLAT):

                for j in range(NLON):

                    for k in range(NTAU - 1):

                        data[i][j][k+1][param_col] = data[i][j][k+1][param_col] + data[i][j][k][param_col]




            # define new half-step grid to interpolate over

            dz = (np.max(z_new) - np.min(z_new)) / (NTAU_new - 1)

            half_grid = np.zeros(NTAU_new + 1)

            half_grid[0] = z_new[0] + (0.5 * dz)

            for n in range(len(half_grid) - 1):

                half_grid[n+1] = half_grid[n] - dz



            # do log interpolation


            for i in range(NLAT):

                for j in range(NLON):


                    # old altitude grid at this lat, lon

                    z_old = data[i][j][:,3]



                    # parameter values on old altitude grid

                    param_old = data[i][j][:,param_col]



                    # make sure no values fall below some epsilon (zeros are bad!)

                    epsilon = 1e-10

                    for n in range(len(param_old)):

                        if param_old[n] < epsilon:

                            param_old[n] = epsilon



                    # log of old parameter

                    log_param_old = np.log(param_old)



                    # linear interpolation function created from old altitudes and log parameter values

                    param_interp = interpolate.interp1d(z_old, log_param_old, kind="linear", fill_value="extrapolate")



                    # apply interpolation function

                    log_param_new = param_interp(half_grid)



                    # exponentiate interpolated valued to get back to original format

                    param_new = np.exp(log_param_new)

                    

                    # discretize back to grid and change parameter values in data array to the new interpolated values

                    for k in range(NTAU_new):

                        data_new[i][j][k][param_col] = param_new[k+1] - param_new[k]





        else:



            for i in range(NLAT):

                for j in range(NLON):


                    # old altitude grid at this lat, lon

                    z_old = data[i][j][:,3]



                    # parameter values on old altitude grid

                    param_old = data[i][j][:,param_col]



                    # make sure no values fall below some epsilon (zeros are bad!)

                    epsilon = 1e-10

                    for n in range(len(param_old)):

                        if param_old[n] < epsilon:

                            param_old[n] = epsilon



                    # log of old parameter

                    log_param_old = np.log(param_old)



                    # linear interpolation function created from old altitudes and log parameter values

                    param_interp = interpolate.interp1d(z_old, log_param_old, kind="linear", bounds_error=False, fill_value=-1e10)



                    # apply interpolation function (values not on new grid should go to zero)

                    log_param_new = param_interp(z_new)



                    # exponentiate interpolated valued to get back to original format

                    param_new = np.exp(log_param_new)



                    # change parameter values in data array to the new interpolated values

                    for k in range(NTAU_new):

                        data_new[i][j][k][param_col] = param_new[k]




        return data_new

        




    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#





    ############################################################
    ###                                                         ###
    ###    ----- MAIN PART OF CODE -- CALL FUNCTIONS HERE ----- ###
    ###                                                         ###
    ############################################################




    # load data and reshape data into more convenient dimensions
    print ('Old file:', old_file)
    data = np.loadtxt(old_file)
    

    # Get the largest value in the 3rd column of the data, pressure
    max_pressure_bar_in_model = data[:,4].max()

    if max_pressure_bar < max_pressure_bar_in_model:
        data = pd.DataFrame(data)
        data = data[data[4] <= max_pressure_bar]
        data = data.to_numpy()
        max_val = int(data[:,2].max())
        NTAU = max_val
        print('cutting off the bottom of the atmosphere')
        
    data = data.reshape((NLAT, NLON, NTAU, NPARAMS))
    data_new = np.zeros((NLAT, NLON, NTAU_new, NPARAMS))


    print(data_new.shape)




    # get new altitude grid

    z_grid = altitudes(data)



    # smooth temperatures (to mitigate numerical noise in cloudy models)

    if smoothing == True:

        for i in range(NLAT):

            for j in range(NLON):

                temps = data[i][j][:,5]

                smooth_temps = savgol_filter(temps, 7, 3)        # (data, window size, polynomial degree)

                #coeffs = np.polyfit(np.arange(len(temps)), temps, 5)  # Fit a 3rd order polynomial
                #interpolated_temps = np.polyval(coeffs, np.arange(len(temps)))  # Evaluate the polynomial
                
                data[i][j][:,5] = smooth_temps



    # interpolate pressures onto new grid (logarithmic)
    data_new = LogInterp_1d(data, data_new, z_grid, 4)


    # interpolate temperature onto new grid (linear)
    data_new = LInterp_1d(data, data_new, z_grid, 5)


    # interpolate winds onto new grid (linear)
    data_new = LInterp_1d(data, data_new, z_grid, 6)
    data_new = LInterp_1d(data, data_new, z_grid, 7)
    data_new = LInterp_1d(data, data_new, z_grid, 8)


    # interpolate optical depths onto new grid (logarithmic)

    if (CLOUDS == 0) and (HAZES == False):
        data_new = LInterp_1d(data, data_new, z_grid, 9) 
        data_new = LInterp_1d(data, data_new, z_grid, 10)
        data_new = LInterp_1d(data, data_new, z_grid, 11)
    else:
        # 1
        data_new = LInterp_1d(data, data_new, z_grid, 9) 
        data_new = LInterp_1d(data, data_new, z_grid, 10)
        data_new = LInterp_1d(data, data_new, z_grid, 11)

        # 2
        data_new = LInterp_1d(data, data_new, z_grid, 12)
        data_new = LInterp_1d(data, data_new, z_grid, 13)
        data_new = LInterp_1d(data, data_new, z_grid, 14)

        # 3
        data_new = LInterp_1d(data, data_new, z_grid, 15)
        data_new = LInterp_1d(data, data_new, z_grid, 16)
        data_new = LInterp_1d(data, data_new, z_grid, 17)

        # 4
        data_new = LInterp_1d(data, data_new, z_grid, 18)
        data_new = LInterp_1d(data, data_new, z_grid, 19)
        data_new = LInterp_1d(data, data_new, z_grid, 20)

        # 5
        data_new = LInterp_1d(data, data_new, z_grid, 21)
        data_new = LInterp_1d(data, data_new, z_grid, 22)
        data_new = LInterp_1d(data, data_new, z_grid, 23)

        # 6
        data_new = LInterp_1d(data, data_new, z_grid, 24)
        data_new = LInterp_1d(data, data_new, z_grid, 25)
        data_new = LInterp_1d(data, data_new, z_grid, 26)

        # 7
        data_new = LInterp_1d(data, data_new, z_grid, 27)
        data_new = LInterp_1d(data, data_new, z_grid, 28)
        data_new = LInterp_1d(data, data_new, z_grid, 29)

        # 8
        data_new = LInterp_1d(data, data_new, z_grid, 30)
        data_new = LInterp_1d(data, data_new, z_grid, 31)
        data_new = LInterp_1d(data, data_new, z_grid, 32)

        # 9
        data_new = LInterp_1d(data, data_new, z_grid, 33)
        data_new = LInterp_1d(data, data_new, z_grid, 34)
        data_new = LInterp_1d(data, data_new, z_grid, 35)

        # 10
        data_new = LInterp_1d(data, data_new, z_grid, 36)
        data_new = LInterp_1d(data, data_new, z_grid, 37)
        data_new = LInterp_1d(data, data_new, z_grid, 38)

        # 11
        data_new = LInterp_1d(data, data_new, z_grid, 39)
        data_new = LInterp_1d(data, data_new, z_grid, 40)
        data_new = LInterp_1d(data, data_new, z_grid, 41)

        # 12
        data_new = LInterp_1d(data, data_new, z_grid, 42)
        data_new = LInterp_1d(data, data_new, z_grid, 43)
        data_new = LInterp_1d(data, data_new, z_grid, 44)

        # 13
        data_new = LInterp_1d(data, data_new, z_grid, 45)
        data_new = LInterp_1d(data, data_new, z_grid, 46)
        data_new = LInterp_1d(data, data_new, z_grid, 47)

        # HAZES
        data_new = LInterp_1d(data, data_new, z_grid, 48)
        data_new = LInterp_1d(data, data_new, z_grid, 49)
        data_new = LInterp_1d(data, data_new, z_grid, 50)


    # lastly, set all altitude grids equal (to new grid) and add lat, lon, level

    for i in range(NLAT):

        for j in range(NLON):

            for k in range(NTAU_new):
                data_new[i][j][k][3] = z_grid[k]
                data_new[i][j][k][2] = k + 1
                data_new[i][j][k][1] = data[i][j][0][1]
                data_new[i][j][k][0] = data[i][j][0][0]


    # double all data, then save to new output file

    np.savetxt(new_file, data_new.reshape(NLAT * NLON * NTAU_new, NPARAMS),
               fmt=' '.join(['%5.4f']*2 + ['%3d']*1 + ['%9.6e']*1 + ['%9.4E']*5 + ['%9.4E']*42 + ['\t']))
    return None



