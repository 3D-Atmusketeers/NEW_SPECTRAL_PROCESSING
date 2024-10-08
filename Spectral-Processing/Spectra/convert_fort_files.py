import numpy as np
import pandas as pd
import os

def convert_to_correct_format(path, runname, planet_name, INITIAL_NTAU, surfp, oom, grav, gasconst, NLAT, NLON):
    def readfortfiles(path,runname,fort26,fort50,nlay,gasconst,oom,grav): 
        with open(path+runname+'/'+fort26) as f:
            first_line=f.readline()
            nlat,nlon,nlev=first_line.split()
            nlat,nlon,nlev=int(nlat),int(nlon),int(nlev)
        f.close()
        data26=np.empty([nlon*nlat*nlev, 6])
        
        l=0
        lp=0
        with open(path+runname+'/'+fort26) as f:
            for line in f:
                if l==0:
                    l+=1
                    continue
                elif l%2==1 and l<=nlon*nlat*nlev*2.:
                    line_pair=np.empty([6])
                    lon, lat, lev, u, v = line.split()
                    line_pair[:5] = np.float32(lon), np.float32(lat), int(lev), np.float32(u), np.float32(v)
                elif l%2==0 and l<=nlon*nlat*nlev*2.:
                    line_pair[5]=np.float32(line)
                    data26[lp,:]=line_pair
                    lp+=1
                elif l>nlon*nlat*nlev*2.:
                    break
                l+=1
        f.close()

        lon_arr_f=data26[:,0]
        lon_arr=np.array([])
        for l in range(0,len(lon_arr_f)):
            el=lon_arr_f[l]
            if not el in lon_arr:
                lon_arr=np.append(lon_arr,el)

        lat_arr_f=data26[:,1]
        lat_arr=np.array([])
        for l in range(0,len(lat_arr_f)):
            el=lat_arr_f[l]
            if not el in lat_arr:
                lat_arr=np.append(lat_arr,el)

        lev_arr_f=data26[:,2]
        lev_arr=np.array([])
        for l in range(0,len(lev_arr_f)):
            el=lev_arr_f[l]
            if not el in lev_arr:
                lev_arr=np.append(lev_arr,el)

        data_26=np.empty([nlev,nlon,nlat,6])
        for l in range(0,data26.shape[0]):
            lon,lat,lev=data26[l,:3]
            lon_i,lat_i,lev_i=np.where(lon_arr==lon)[0][0],np.where(lat_arr==lat)[0][0],np.where(lev_arr==lev)[0][0]
            data_26[lev_i,lon_i,lat_i,:]=data26[l,:]

        nlev,nlon,nlat,nparam=data_26.shape
        temps=data_26[:,:,:,5]
        
        #get surface pressures
        with open(path+runname+'/'+fort50,'r') as data_50:  #fort_50 long1 lat1 pressure 
            specificp=np.zeros((NLAT,NLON))  #lat,long               #long1 lat2 pressure
            acount=0
            bcount=0
            next(data_50)
            for line in data_50:
                p=line.split()
                if bcount<NLAT and acount<NLON:
                    specificp[bcount][acount]=((float(p[2]))+1)*surfp
                    bcount=bcount+1
                else:
                    if acount<NLON:
                        acount=acount+1
                        bcount=0
                        if acount<NLON:
                            specificp[bcount][acount]=((float(p[2]))+1)*surfp
                            bcount=bcount+1

        sigma=np.empty([nlay])*0.0
        if oom>0: #setting up pressure values 
            stp=-1.0*oom/nlay
            sigma[nlay-1]=10.**(stp/2.)
            for n in range(nlay-2,-1,-1):
                sigma[n]=sigma[n+1]*10.**(stp)

        p_BAR=sigma*surfp
        sp=specificp

        z=np.zeros((nlat,nlon,nlay,2))#one for z value, one for prssure
        z[:,:,nlay-1,0]=(gasconst/grav) *.5*(temps[nlay-1,:,:].T) * np.log(surfp/sp/sigma[nlay-1])
        z[:,:,-1,1]=p_BAR[-1]


        #integrate hydrostatic to solve for higher levels
        start=nlay-2
        while start >= 0: #matches idl
            z[:,:,start,0]=z[:,:,start+1,0] + (gasconst/grav) *0.5*(temps[start,:,:].T+temps[start+1,:,:].T) * np.log(sigma[start+1]/sigma[start])
            z[:,:,start,1]=p_BAR[start] 
            start=start - 1

        return data_26,nlon,nlat,nlev,nparam,z


    levs=INITIAL_NTAU
    if os.path.isfile(path+runname+'/fort.2600'):
        data_26,nlon,nlat,nlev,nparam,z=readfortfiles(path,runname,'fort.2600','fort.5000',levs,gasconst,oom,grav)
        print ("Using the fort.2600 and fort.5000 files")
    else:
        data_26,nlon,nlat,nlev,nparam,z=readfortfiles(path,runname,'fort.26','fort.50',levs,gasconst,oom,grav)
        print ("Using the fort.26 and fort.50 files")


    # Reshape the array to be 2D
    # It looks like there are INITIAL_NTAU levels, 48 lats and 96 lons
    num_variables = 6
    df = data_26.reshape(levs * NLON * NLAT, num_variables)

    # Make it a pandas dataframe
    pd_df = pd.DataFrame(df, columns=['lon', 'lat', 'level', 'u', 'v', 'temps'])

    # Make the z stuff also a dataframe
    z_df = pd.DataFrame(z.reshape(levs * NLON * NLAT, 2), columns=['alt','pressure'])

    # Sort all the values by lat, then lon, then the level
    data = pd_df.sort_values(by=['lat', 'lon', 'level'], axis=0, ascending=[True, True, True])

    # Reset the file index and drop the resulting index columns
    data = data.reset_index(drop=True)

    # Concatenate the two
    data = pd.concat([data, z_df], axis=1)

    #winds speeds at boundaries are weird
    data['u'][data['level'] == 1] = 0
    data['v'][data['level'] == 1] = 0

    data['w'] = 0

    # Reshuffle the order of the lats and the lons
    data = data[['lat', 'lon', 'level','alt','pressure', 'temps', 'u', 'v', 'w']]

    data = data.sort_values(by=['lat', 'lon', 'level'], axis=0, ascending=[False, True, True])


    # Save the reformatted data
    np.savetxt('../PLANET_MODELS/' + planet_name + '.txt', data.values, fmt='%5.4f %6.4f %3d %9.4E %9.4E %9.4E %9.4E %9.4E %9.4E')


