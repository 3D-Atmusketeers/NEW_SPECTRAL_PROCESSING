import pandas as pd
import numpy as np
import os
from tqdm import tqdm
pd.options.mode.chained_assignment = None  # default='warn'


def add_clouds_to_gcm_output(path, runname, planet_name, grav, MTLX, CLOUDS, MOLEF,
                             aerosol_layers, INITIAL_NTAU, GASCON, HAZE_TYPE, HAZES, wav_loc, MET_X_SOLAR, path_to_data='DATA/'):
    column_names = ['lat' , 'lon', 'level' , 'altitude(m)',
                    'pressure(bars)', 'temp(k)',
                    'EW vel(m/s)','NS vel','vert vel']
    df = pd.read_csv('../PLANET_MODELS/' + planet_name + '.txt', delim_whitespace=True, names=column_names, index_col=False)

    df['tau1'] = 0.
    df['g01']  = 0.
    df['pi01'] = 0.

    df['tau2'] = 0.
    df['g02']  = 0.
    df['pi02'] = 0.

    df['tau3'] = 0.
    df['g03']  = 0.
    df['pi03'] = 0.

    df['tau4'] = 0.
    df['g04']  = 0.
    df['pi04'] = 0.

    df['tau5'] = 0.
    df['g05']  = 0.
    df['pi05'] = 0.
    
    df['tau6'] = 0.
    df['g06']  = 0.
    df['pi06'] = 0.
    
    df['tau7'] = 0.
    df['g07']  = 0.
    df['pi07'] = 0.
    
    df['tau8'] = 0.
    df['g08']  = 0.
    df['pi08'] = 0.
    
    df['tau9'] = 0.
    df['g09']  = 0.
    df['pi09'] = 0.
    
    df['tau10'] = 0.
    df['g010']  = 0.
    df['pi010'] = 0.
    
    df['tau11'] = 0.
    df['g011']  = 0.
    df['pi011'] = 0.
    
    df['tau12'] = 0.
    df['g012']  = 0.
    df['pi012'] = 0.
    
    df['tau13'] = 0.
    df['g013']  = 0.
    df['pi013'] = 0.

    df['tau_haze'] = 0.0
    df['g0_haze']  = 0.0
    df['pi0_haze'] = 0.0

    df_copy = df.copy(deep=True)

    # The input particle size is used ONLY for the lookup of the scattering properties in the scattering tables
    input_pressure_array_cgs               = np.asarray([1.000e+00, 1.460e+00, 2.120e+00, 3.090e+00, 4.500e+00, 6.550e+00, 9.540e+00, 1.389e+01, 2.024e+01, 2.947e+01, 4.292e+01, 6.251e+01, 9.103e+01, 1.326e+02, 1.931e+02, 2.812e+02, 4.095e+02, 5.964e+02, 8.685e+02, 1.265e+03, 1.842e+03, 2.683e+03, 3.907e+03, 5.690e+03, 8.286e+03, 1.207e+04, 1.758e+04, 2.560e+04, 3.728e+04, 5.429e+04, 7.906e+04, 1.151e+05, 1.677e+05, 2.442e+05, 3.556e+05, 5.179e+05, 7.543e+05, 1.099e+06, 1.600e+06, 2.330e+06, 3.393e+06, 4.942e+06, 7.197e+06, 1.048e+07, 1.526e+07, 2.223e+07, 3.237e+07, 4.715e+07, 6.866e+07, 1.000e+08])
    input_particle_size_array_in_meters    = np.asarray([1.000e-08, 1.207e-08, 1.456e-08, 1.758e-08, 2.121e-08, 2.560e-08, 3.089e-08, 3.728e-08, 4.498e-08, 5.429e-08, 6.551e-08, 7.906e-08, 9.541e-08, 1.151e-07, 1.389e-07, 1.677e-07, 2.024e-07, 2.442e-07, 2.947e-07, 3.556e-07, 4.292e-07, 5.179e-07, 6.251e-07, 7.543e-07, 9.103e-07, 1.099e-06, 1.326e-06, 1.600e-06, 1.931e-06, 2.330e-06, 2.812e-06, 3.393e-06, 4.095e-06, 4.942e-06, 5.964e-06, 7.197e-06, 8.685e-06, 1.048e-05, 1.265e-05, 1.526e-05, 1.842e-05, 2.223e-05, 2.683e-05, 3.237e-05, 3.907e-05, 4.715e-05, 5.690e-05, 6.866e-05, 8.286e-05, 1.000e-04])
    particle_size_vs_layer_array_in_meters = np.asarray([1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.000e-07, 1.007e-07, 1.048e-07, 1.111e-07, 1.201e-07, 1.333e-07, 1.556e-07, 2.106e-07, 2.655e-07, 3.452e-07, 4.609e-07, 6.302e-07, 8.764e-07, 1.235e-06, 1.757e-06, 2.517e-06, 3.624e-06, 5.237e-06, 7.585e-06, 1.100e-05, 1.599e-05, 2.324e-05, 3.380e-05, 4.918e-05, 7.159e-05, 1.042e-04])

    DENSITY = [1.98e3,4.09e3,1.86e3,4.0e3,5.22e3,2.65e3,3.27e3,5.76e3,8.9e3,7.9e3,3.34e3,3.98e3,3.95e3]

    CLOUD_MOLAR_MASSES = [74.55E-3,   #
                          97.47E-3,   # ZnS
                          78.05E-3,   # Na2S
                          87.00E-3,   # MnS
                          52.00E-3,   # Cr
                          60.08E-3,   # SiO2
                          160.95E-3,  # Mg2Si04
                          66.94E-3,   # VO
                          58.69E-3,   # Ni
                          55.85E-3,   # Fe
                          172.23E-3,  # Ca2Si04
                          135.94E-3,  # CaTiO3
                          102.00E-3]  # Al2O3

    haze_pressure_array_pascals = [9.999999999999999167e-02,1.204503540258781008e-01,1.450828778495940330e-01,
                                   1.747528400007682947e-01,2.104904144512021735e-01,2.535364493970111432e-01,
                                   3.053855508833412391e-01,3.678379771828634293e-01,4.430621457583877598e-01,
                                   5.336699231206313288e-01,6.428073117284319737e-01,7.742636826811276629e-01,
                                   9.326033468832199969e-01,1.123324032978026521e+00,1.353047774579807516e+00,
                                   1.629750834620643518e+00,1.963040650040272617e+00,2.364489412645407018e+00,
                                   2.848035868435804918e+00,3.430469286314919319e+00,4.132012400115334216e+00,
                                   4.977023564332114347e+00,5.994842503189408589e+00,7.220809018385470957e+00,
                                   8.697490026177835176e+00,1.047615752789665144e+01,1.261856883066021062e+01,
                                   1.519911082952934755e+01,1.830738280295369691e+01,2.205130739903045622e+01,
                                   2.656087782946686460e+01,3.199267137797384564e+01,3.853528593710531425e+01,
                                   4.641588833612782139e+01,5.590810182512228010e+01,6.734150657750828373e+01,
                                   8.111308307896872805e+01,9.770099572992256753e+01,1.176811952434999142e+02,
                                   1.417474162926806400e+02,1.707352647470692091e+02,2.056512308348653448e+02,
                                   2.477076355991711409e+02,2.983647240283340238e+02,3.593813663804629073e+02,
                                   4.328761281083061476e+02,5.214008287999689628e+02,6.280291441834259558e+02,
                                   7.564633275546291316e+02,9.111627561154896284e+02,1.097498765493056908e+03,
                                   1.321941148466031336e+03,1.592282793341094020e+03,1.917910261672488559e+03,
                                   2.310129700083162788e+03,2.782559402207125913e+03,3.351602650938848001e+03,
                                   4.037017258596557895e+03,4.862601580065353119e+03,5.857020818056673306e+03,
                                   7.054802310718645458e+03,8.497534359086455879e+03,1.023531021899026928e+04,
                                   1.232846739442068429e+04,1.484968262254466572e+04,1.788649529057434847e+04,
                                   2.154434690031886566e+04,2.595024211399737396e+04,3.125715849688241360e+04,
                                   3.764935806792471703e+04,4.534878508128591056e+04,5.462277217684347852e+04,
                                   6.579332246575682075e+04,7.924828983539185720e+04,9.545484566618347890e+04,
                                   1.149756995397738065e+05,1.384886371393874579e+05,1.668100537200059334e+05,
                                   2.009233002565049974e+05,2.420128264794383431e+05,2.915053062825181987e+05,
                                   3.511191734215134638e+05,4.229242874389507924e+05,5.094138014816385694e+05,
                                   6.135907273413175717e+05,7.390722033525790321e+05,8.902150854450392071e+05,
                                   1.072267222010325408e+06,1.291549665014885366e+06,1.555676143930475460e+06,
                                   1.873817422860386781e+06,2.257019719633921515e+06,2.718588242732945830e+06,
                                   3.274549162877731491e+06,3.944206059437664226e+06,4.750810162102803588e+06,
                                   5.722367659350220114e+06,6.892612104349709116e+06,8.302175681319752708e+06,
                                   1.000000000000000000e+07]

    FMOLW = [0] * len(CLOUD_MOLAR_MASSES)
    GAS_CONSTANT_R = 8.314462618
    for j in range(len(CLOUD_MOLAR_MASSES)):
        FMOLW[j] = CLOUD_MOLAR_MASSES[j] / (GAS_CONSTANT_R/GASCON)

    CORFACT = [0.005,0.018,0.050,0.135,.367,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000]
    
    if (0.1  <= MET_X_SOLAR < 5.0):
        metstr = '1xsolar'
    elif (5.0 <= MET_X_SOLAR < 20.0):
        metstr = '10xsolar'
    elif (20.0 <= MET_X_SOLAR < 50.0):
        metstr = '30xsolar'
    elif (50.0 <= MET_X_SOLAR < 150.0):
        metstr = '100xsolar'
    elif (150.0 <= MET_X_SOLAR):
        metstr = '300xsolar'
    else:
        print("Error in setting condensation curves")
    Tconds = np.loadtxt(os.path.join('DATA', 'condcurves', f'condcurves_{metstr}_new.txt'), usecols=range(1,14)).T # first column is pressure axis

    qe0  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'KCl_wav_qext.txt'))
    g00  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'KCl_wav_gg.txt'))
    pi00 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'KCl_wav_pi0.txt'))

    qe1  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'ZnS_wav_qext.txt'))
    g01  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'ZnS_wav_gg.txt'))
    pi01 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'ZnS_wav_pi0.txt'))

    qe2  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Na2S_wav_qext.txt'))
    g02  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Na2S_wav_gg.txt'))
    pi02 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Na2S_wav_pi0.txt'))

    qe3  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'MnS_wav_qext.txt'))
    g03  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'MnS_wav_gg.txt'))
    pi03 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'MnS_wav_pi0.txt'))

    qe4  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Cr_wav_qext.txt'))
    g04  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Cr_wav_gg.txt'))
    pi04 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Cr_wav_pi0.txt'))

    qe5  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'SiO2_wav_qext.txt'))
    g05  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'SiO2_wav_gg.txt'))
    pi05 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'SiO2_wav_pi0.txt'))

    qe6  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Mg2SiO4_wav_qext.txt'))
    g06  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Mg2SiO4_wav_gg.txt'))
    pi06 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Mg2SiO4_wav_pi0.txt'))

    qe7  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'VO_wav_qext.txt'))
    g07  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'VO_wav_gg.txt'))
    pi07 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'VO_wav_pi0.txt'))

    qe8  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Ni_wav_qext.txt'))
    g08  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Ni_wav_gg.txt'))
    pi08 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Ni_wav_pi0.txt'))

    qe9  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Fe_wav_qext.txt'))
    g09  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Fe_wav_gg.txt'))
    pi09 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Fe_wav_pi0.txt'))

    qe10  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaSiO4_wav_qext.txt'))
    g010  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaSiO4_wav_gg.txt'))
    pi010 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaSiO4_wav_pi0.txt'))

    qe11  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaTiO3_wav_qext.txt'))
    g011  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaTiO3_wav_gg.txt'))
    pi011 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'CaTiO3_wav_pi0.txt'))

    qe12  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Al2O3_wav_qext.txt'))
    g012  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Al2O3_wav_gg.txt'))
    pi012 = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'Al2O3_wav_pi0.txt'))

    if (HAZES == True):
        print ('haze_' + HAZE_TYPE + '_wav_tau_per_bar.txt')
        tau_haze = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_' + HAZE_TYPE + '_wav_tauperbar.txt'))
        g0_haze  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_'  + HAZE_TYPE  + '_wav_gg.txt'))
        pi0_haze = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_'  + HAZE_TYPE  + '_wav_pi0.txt'))
    else:
        tau_haze = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_tholin_wav_tauperbar.txt'))
        g0_haze  = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_tholin_wav_gg.txt'))
        pi0_haze = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'haze_tholin_wav_pi0.txt'))

        tau_haze = tau_haze * 0
        g0_haze  = g0_haze  * 0
        pi0_haze = pi0_haze * 0
        print ('No hazes')
        pass


    QE_OPPR  = np.array([qe0, qe1,qe2,qe3,qe4,qe5,qe6,qe7,qe8,qe9,qe10,qe11,qe12,tau_haze], dtype=object)
    PI0_OPPR = np.array([pi00, pi01,pi02,pi03,pi04,pi05,pi06,pi07,pi08,pi09,pi010,pi011,pi012,pi0_haze], dtype=object)
    G0_OPPR  = np.array([g00, g01,g02,g03,g04,g05,g06,g07,g08,g09,g010,g011,g012,g0_haze], dtype=object)
    
    G = grav

    pressure_array_for_cloud_scattering_data_in_pascals = np.loadtxt(os.path.join(path_to_data, 'Aerosol_Data', 'pressure_array_for_cloud_scattering_data_in_pascals.txt'))

    print ("Adding Clouds, no scattering params, this just fills in 0s unless MOLEF is specified")
    
    ### EVERY LINE ABOVE THIS ONE HAS NOT BEEN EDITED BY THOMAS
    
    nlat, nlon, nlay = len(df.lat.unique()),len(df.lon.unique()),INITIAL_NTAU
    
    ### Defining all the same constants as in old version but as nlay-shaped arrays
    pres = df.sort_values(by=['lat','lon','pressure(bars)'])['pressure(bars)'].values.reshape(nlat,nlon,nlay)[0][0]
    layer_index = np.array([np.abs(input_pressure_array_cgs - pres[i]*1e6).argmin() for i in range(len(pres))])
    particle_size = particle_size_vs_layer_array_in_meters[layer_index]
    
    delta_pres = np.array([(pres[i+1]-pres[i])/2 if i==0 else (pres[i]-pres[i-1]) for i in range(len(pres))])
    dpg = 1e5/grav*delta_pres
    
    temp = df.sort_values(by=['lat','lon','pressure(bars)'])['temp(k)'].values.reshape(nlat,nlon,nlay)
    size_loc = [np.abs(input_particle_size_array_in_meters  - particle_size[i]).argmin() for i in range(len(particle_size))]

    haze_layer_index = [np.abs(haze_pressure_array_pascals - pres[i]*1e5).argmin() for i in range(len(pres))]
    
    CORFACT = np.array(CORFACT)[layer_index]
    pressure_for_scattering_data_loc = [np.abs(pressure_array_for_cloud_scattering_data_in_pascals - pres[i]*1e5).argmin() for i in range(len(pres))]
    df = df.sort_values(by=['lat','lon','pressure(bars)'])
    
    for cidx in range(13): # loop through species
        tcond = Tconds[cidx][layer_index] # cond Ts for each layer
        condfact = temp[:,:] <= tcond # everywhere below the cond curve

#       condfact[:,:,0] = False # delete clouds from top layer
        botlay = condfact.cumsum(axis=2).argmax(axis=2).ravel() # find index of cloud base
        
        ### mask out everything above the cloud top
        laymaxmsk = np.array([],dtype=int)
        for i in range(int(nlat*nlon)):
            laymaxmsk = np.append(laymaxmsk, i*nlay + np.array(range(int(botlay[i]+1-aerosol_layers))))
        np.put(condfact,laymaxmsk.astype(int),0.0)

        ### Do the "if close to condensation curve, fraction condenses" calculation
        partcld = np.where((condfact.ravel()) & (temp[:,:] > tcond-10).ravel())[0]
        condfact = condfact.astype(float) # all the places with clouds are now 1, all others are 0
        np.put(condfact,partcld,((tcond - temp[:,:])/10).ravel()[partcld])
        
        ### Calculate the actual cloud properties
        tau = dpg * MOLEF[cidx]*3./4./particle_size/DENSITY[cidx]*FMOLW[cidx]*condfact*MTLX*CORFACT
        tau_wl_dep = tau * QE_OPPR[cidx][pressure_for_scattering_data_loc][:,wav_loc]
        
        ### Everywhere clouds exist (even partially), set scattering params
        tau2 = tau.copy()
        tau2[tau>0] = 1.0
        g0 = tau2 * G0_OPPR[cidx][size_loc][:,wav_loc]
        pi0 = tau2 * PI0_OPPR[cidx][size_loc][:,wav_loc]
        
        ### Save to dfs
        df['tau'+str(cidx+1)] = tau.ravel()
        df_copy['tau'+str(cidx+1)] = tau_wl_dep.ravel()
        df_copy['g0'+str(cidx+1)] = g0.ravel()
        df_copy['pi0'+str(cidx+1)] = pi0.ravel()
        

    if HAZES == True:
        df['tau_haze'] = 1.0
        df['g0_haze']  = 1.0
        df['pi0_haze'] = 1.0

        df_copy['tau_haze'] = (np.ones(np.shape(tau)) * delta_pres * QE_OPPR[13][haze_layer_index][:,wav_loc]).ravel()
        df_copy['g0_haze']  = (np.ones(np.shape(tau)) * G0_OPPR[13][haze_layer_index][:,wav_loc]).ravel()
        df_copy['pi0_haze'] = (np.ones(np.shape(tau)) * PI0_OPPR[13][haze_layer_index][:,wav_loc]).ravel()

        
    planet_file_with_clouds = '../PLANET_MODELS/' + planet_name + '_with_clouds.txt'
    np.savetxt(planet_file_with_clouds, df.values, fmt=' '.join(['%5.4f']*2 + ['%3d']*1 + ['%9.6E']*1 + ['%9.4E']*5 + ['%9.4E']*42 + ['\t']))
    
    planet_file_with_clouds = '../PLANET_MODELS/' + planet_name + '_with_clouds_and_wavelength_dependence.txt'
    np.savetxt(planet_file_with_clouds, df_copy.values, fmt=' '.join(['%5.4f']*2 + ['%3d']*1 + ['%9.6E']*1 + ['%9.4E']*5 + ['%9.4E']*42 + ['\t']))
