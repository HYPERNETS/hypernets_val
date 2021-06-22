#!/usr/bin/env python3
# coding: utf-8
"""
Created on Tue Jun 22 15:09:33 2021

@author: javier.concha
"""
"""
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
def create_extract(size_box,station_name,path_source,path_output,in_situ_lat,in_situ_lon,res_str,insitu_sensor):
    # extract IFP-OL-2 version
    with open(os.path.join(path_source,'xfdumanifest.xml'),'r', encoding="utf-8") as read_obj:
        check_version = False
        for line in read_obj:
            if 'IPF-OL-2' in line and check_version == False:
                IPF_OL_2_version = line.split('"')[3]
                proc_version_str = f'IPF-OL-2 version {IPF_OL_2_version}'
                print(proc_version_str)
                check_version = True
                pass

    
    #% open nc file
    coordinates_filename = 'geo_coordinates.nc'

    rhow_0400p00_filename = 'Oa01_reflectance.nc'
    rhow_0412p50_filename = 'Oa02_reflectance.nc'
    rhow_0442p50_filename = 'Oa03_reflectance.nc'
    rhow_0490p00_filename = 'Oa04_reflectance.nc'
    rhow_0510p00_filename = 'Oa05_reflectance.nc'
    rhow_0560p00_filename = 'Oa06_reflectance.nc'
    rhow_0620p00_filename = 'Oa07_reflectance.nc'
    rhow_0665p00_filename = 'Oa08_reflectance.nc'
    rhow_0673p75_filename = 'Oa09_reflectance.nc'
    rhow_0681p25_filename = 'Oa10_reflectance.nc'
    rhow_0708p75_filename = 'Oa11_reflectance.nc'
    rhow_0753p75_filename = 'Oa12_reflectance.nc'
    rhow_0778p75_filename = 'Oa16_reflectance.nc'
    rhow_0865p00_filename = 'Oa17_reflectance.nc'
    rhow_0885p00_filename = 'Oa18_reflectance.nc'
    rhow_1020p50_filename = 'Oa21_reflectance.nc'
    
    AOT_0865p50_filename = 'w_aer.nc'
    WQSF_filename = 'wqsf.nc'
    
    filepah = os.path.join(path_source,coordinates_filename)
    nc_sat = Dataset(filepah,'r')
    
    lat = nc_sat.variables['latitude'][:,:]
    lon = nc_sat.variables['longitude'][:,:]
    
    contain_flag = cfs.contain_location(lat,lon,in_situ_lat,in_situ_lon)
    
    if contain_flag:
        if not args.verbose:
            print('-----------------')
        r, c = cfs.find_row_column_from_lat_lon(lat,lon,in_situ_lat,in_situ_lon)
        
        start_idx_x = (r-int(size_box/2))
        stop_idx_x = (r+int(size_box/2)+1)
        start_idx_y = (c-int(size_box/2))
        stop_idx_y = (c+int(size_box/2)+1)
    
        
        if r>=0 and r+1<lat.shape[0] and c>=0 and c+1<lat.shape[1]:
            # read nc file
            filepah = os.path.join(path_source,rhow_0400p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0400p00 = nc_sat.variables['Oa01_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source,rhow_0412p50_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0412p50 = nc_sat.variables['Oa02_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0442p50_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0442p50 = nc_sat.variables['Oa03_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0490p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0490p00 = nc_sat.variables['Oa04_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0510p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0510p00 = nc_sat.variables['Oa05_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0560p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0560p00 = nc_sat.variables['Oa06_reflectance'][:]
            nc_sat.close()
    
            filepah = os.path.join(path_source,rhow_0620p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0620p00 = nc_sat.variables['Oa07_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0665p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0665p00 = nc_sat.variables['Oa08_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_0673p75_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0673p75 = nc_sat.variables['Oa09_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source,rhow_0681p25_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0681p25 = nc_sat.variables['Oa10_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source,rhow_0708p75_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0708p75 = nc_sat.variables['Oa11_reflectance'][:]
            nc_sat.close() 

            filepah = os.path.join(path_source,rhow_0753p75_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0753p75 = nc_sat.variables['Oa12_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source,rhow_0778p75_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0778p75 = nc_sat.variables['Oa16_reflectance'][:]
            nc_sat.close()          
            
            filepah = os.path.join(path_source,rhow_0865p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0865p00 = nc_sat.variables['Oa17_reflectance'][:]
            nc_sat.close()

            filepah = os.path.join(path_source,rhow_0885p00_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_0885p00 = nc_sat.variables['Oa18_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,rhow_1020p50_filename)
            nc_sat = Dataset(filepah,'r')
            rhow_1020p50 = nc_sat.variables['Oa21_reflectance'][:]
            nc_sat.close()
            
            filepah = os.path.join(path_source,AOT_0865p50_filename)
            nc_sat = Dataset(filepah,'r')
            AOT_0865p50 = nc_sat.variables['T865'][:]
            nc_sat.close()
    
            filepah = os.path.join(path_source,WQSF_filename)
            nc_sat = Dataset(filepah,'r')
            WQSF = nc_sat.variables['WQSF'][:]
            WQSF_flag_masks = nc_sat.variables['WQSF'].flag_masks
            WQSF_flag_meanings = nc_sat.variables['WQSF'].flag_meanings
            nc_sat.close()

            #%% Calculate BRDF
            ws0, ws1, sza, saa, vza, vaa = extract_wind_and_angles(path_source,in_situ_lat,in_situ_lon)
            
            filepah = os.path.join(path_source,'chl_oc4me.nc')
            nc_sat1 = Dataset(filepah,'r')
            CHL_OC4ME = nc_sat1.variables['CHL_OC4ME'][:]
            nc_sat1.close()
            CHL_OC4ME_extract = ma.array(CHL_OC4ME[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])
            BRDF0 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF1 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF2 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF3 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF4 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF5 = np.full(CHL_OC4ME_extract.shape,np.nan)
            BRDF6 = np.full(CHL_OC4ME_extract.shape,np.nan)
            for ind0 in range(CHL_OC4ME_extract.shape[0]):
                for ind1 in range(CHL_OC4ME_extract.shape[1]):
                    chl = CHL_OC4ME_extract[ind0,ind1]
                    # 412.5, 442.5, 490, 510, 560, 620, 660 bands
                    # 0      1      2    3    4    5    6   brdf index
                    # 412.5  442.5  490  510  560  620  665 OLCI bands
                    # 02     03     04   05   06   07   08  OLCI band names in L2
                    brdf_coeffs = brdf.brdf(ws0, ws1, chl, sza, saa, vza, vaa)
                    BRDF0[ind0,ind1] = brdf_coeffs[0,0]
                    BRDF1[ind0,ind1] = brdf_coeffs[0,1]
                    BRDF2[ind0,ind1] = brdf_coeffs[0,2]
                    BRDF3[ind0,ind1] = brdf_coeffs[0,3]
                    BRDF4[ind0,ind1] = brdf_coeffs[0,4]
                    BRDF5[ind0,ind1] = brdf_coeffs[0,5]
                    BRDF6[ind0,ind1] = brdf_coeffs[0,6]

            #%% Save extract as netCDF4 file
            path_out = os.path.join(path_output,'EXTRACTS')
            filename = path_source.split('/')[-1].replace('.','_')+'_extract_'+station_name+'.nc'
            ofname = os.path.join(path_out,filename)
            
            print(filename)

            satellite = filename[0:2]
            platform = filename[2]
            sensor = 'olci'
            
            if os.path.exists(ofname):
              os.remove(ofname)
            
            new_MDB = Dataset(ofname, 'w', format='NETCDF4')
            new_MDB.MDB_software_version = '0.0'
            new_MDB.creation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")            
            new_MDB.satellite = satellite
            new_MDB.platform = platform
            new_MDB.sensor = sensor
            new_MDB.description = f'{satellite}{platform} {sensor.upper()} {res_str} L2 - {insitu_sensor} Matchup Data Base'
            # new_MDB.satellite_start_time = nc_sat.start_time
            # new_MDB.satellite_stop_time = nc_sat.stop_time    
            # new_MDB.satellite_PDU = path_source.split('/')[-1]
            # new_MDB.satellite_path_source = path_source
            new_MDB.satellite_aco_processor = 'Atmospheric Correction processor: xxx'
            new_MDB.satellite_proc_version = proc_version_str

            new_MDB.datapolicy = 'Notice to users: Add data policy'
            new_MDB.insitu_sensor_processor_version = '0.0'
            new_MDB.insitu_site_name = station_name

            new_MDB.insitu_lat = in_situ_lat
            new_MDB.insitu_lon = in_situ_lon

            new_MDB.satellite_ws0 = ws0
            new_MDB.satellite_ws1 = ws1
            new_MDB.satellite_SZA_center_pixel = sza
            new_MDB.satellite_SAA_center_pixel = saa
            new_MDB.satellite_VZA_center_pixel = vza
            new_MDB.satellite_VAA_center_pixel = vaa
            
            # dimensions
            new_MDB.createDimension('satellite_id', None)
            new_MDB.createDimension('rows', size_box)
            new_MDB.createDimension('columns', size_box)
            new_MDB.createDimension('satellite_bands', 16)
            new_MDB.createDimension('satellite_BRDF_Bands', 7)
            
            
            # variables  
            # satellite_SZA = new_MDB.createVariable('satellite_SZA', 'f4', ('rows','columns'), fill_value=-999, zlib=True, complevel=6)
            # satellite_SZA[:] = [SZA[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
            # satellite_SZA.long_name = 'Sun Zenith Angle'
            # satellite_SZA.long_name = 'Sun Zenith Angle'    
            # satellite_SZA.units = 'degrees'
# SZA
# SAA
# OZA
# OAA

            satellite_time = new_MDB.createVariable('satellite_time',  'f4', ('satellite_id'), fill_value=-999, zlib=True, complevel=6)  
            satellite_time[0] = float(datetime.strptime(nc_sat.start_time,"%Y-%m-%dT%H:%M:%S.%fZ").timestamp())
            satellite_time.units = "Seconds since 1970-1-1"

            satellite_PDU = new_MDB.createVariable('satellite_PDU',  'S2', ('satellite_id'), zlib=True, complevel=6) # string
            satellite_PDU[0] = path_source.split('/')[-1]
            satellite_PDU.long_name = "OLCI source PDU name"

            satellite_latitude = new_MDB.createVariable('satellite_latitude',  'f4', ('satellite_id','rows','columns'), fill_value=-999, zlib=True, complevel=6) 
            satellite_latitude[0,:,:] = [lat[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]
            
            satellite_longitude = new_MDB.createVariable('satellite_longitude',  'f4', ('satellite_id','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_longitude[0,:,:] = [lon[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]]

            # double satellite_bands          (satellite_bands) ;
            satellite_bands = new_MDB.createVariable('satellite_bands',  'f4', ('satellite_bands'), fill_value=-999, zlib=True, complevel=6) 
            satellite_bands[:] = [0400.00,0412.50,0442.50,0490.00,0510.00,0560.00,0620.00,0665.00,0673.75,0681.25,0708.75,0753.75,0778.75,0865.00,0885.00,1020.50]

            # double satellite_BRDF_Bands     (satellite_BRDF_Bands) ;
            satellite_BRDF_Bands = new_MDB.createVariable('satellite_BRDF_Bands',  'f4', ('satellite_BRDF_Bands'), fill_value=-999, zlib=True, complevel=6) 
            satellite_BRDF_Bands[:] = [412.50,442.50,490.00,510.00,560.00,620.00,665.00]
    
            # NOT BRDF-corrected
            satellite_rhow = new_MDB.createVariable('satellite_rhow', 'f4', ('satellite_id','satellite_bands','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_rhow[0,0,:,:] = [ma.array(rhow_0400p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,1,:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,2,:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,3,:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,4,:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,5,:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,6,:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,7,:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,8,:,:] = [ma.array(rhow_0673p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,9,:,:] = [ma.array(rhow_0681p25[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,10,:,:] = [ma.array(rhow_0708p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,11,:,:] = [ma.array(rhow_0753p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,12,:,:] = [ma.array(rhow_0778p75[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,13,:,:] = [ma.array(rhow_0865p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,14,:,:] = [ma.array(rhow_0885p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow[0,15,:,:] = [ma.array(rhow_1020p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_rhow.description = 'Satellite rhow.'
            
            # BRDF-corrected
            satellite_BRDF_rhow = new_MDB.createVariable('satellite_BRDF_rhow', 'f4', ('satellite_id','satellite_BRDF_Bands','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_BRDF_rhow[0,0,:,:] = [ma.array(rhow_0412p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF0)]
            satellite_BRDF_rhow[0,1,:,:] = [ma.array(rhow_0442p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF1)]
            satellite_BRDF_rhow[0,2,:,:] = [ma.array(rhow_0490p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF2)]
            satellite_BRDF_rhow[0,3,:,:] = [ma.array(rhow_0510p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF3)]
            satellite_BRDF_rhow[0,4,:,:] = [ma.array(rhow_0560p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF4)]
            satellite_BRDF_rhow[0,5,:,:] = [ma.array(rhow_0620p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF5)]
            satellite_BRDF_rhow[0,6,:,:] = [ma.array(rhow_0665p00[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y]*BRDF6)]
            satellite_BRDF_rhow.description = 'Satellite rhow BRDF-corrected'
            
            satellite_AOT_0865p50_box = new_MDB.createVariable('satellite_AOT_0865p50', 'f4', ('satellite_id','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_AOT_0865p50_box[0,:,:] = [ma.array(AOT_0865p50[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'
    
            satellite_WQSF = new_MDB.createVariable('satellite_WQSF', 'f4', ('satellite_id','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_WQSF[0,:,:] = [ma.array(WQSF[start_idx_x:stop_idx_x,start_idx_y:stop_idx_y])]
            satellite_WQSF.description = 'Satellite Level 2 WATER Product, Classification, Quality and Science Flags Data Set'
            satellite_WQSF.flag_masks = WQSF_flag_masks
            satellite_WQSF.flag_meanings = WQSF_flag_meanings
            
            satellite_BRDF_fQ = new_MDB.createVariable('satellite_BRDF_fQ', 'f4', ('satellite_id','satellite_BRDF_Bands','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_BRDF_fQ[0,0,:,:] = [ma.array(BRDF0)]
            satellite_BRDF_fQ[0,1,:,:] = [ma.array(BRDF1)]
            satellite_BRDF_fQ[0,2,:,:] = [ma.array(BRDF2)]
            satellite_BRDF_fQ[0,3,:,:] = [ma.array(BRDF3)]
            satellite_BRDF_fQ[0,4,:,:] = [ma.array(BRDF4)]
            satellite_BRDF_fQ[0,5,:,:] = [ma.array(BRDF5)]
            satellite_BRDF_fQ[0,6,:,:] = [ma.array(BRDF6)]
            satellite_BRDF_fQ.description = 'Satellite BRDF fQ coefficients'

            satellite_chl_oc4me = new_MDB.createVariable('chl_oc4me', 'f4', ('satellite_id','rows','columns'), fill_value=-999, zlib=True, complevel=6)
            satellite_chl_oc4me[0,:,:] = [ma.array(CHL_OC4ME_extract)]
            satellite_chl_oc4me.description = 'Satellite Chlorophyll-a concentration from OC4ME.'

            new_MDB.close()
            # print('Extract created!')
                
        else:
            print('Index out of bound!')
    else:
        if args.verbose:
            print('File does NOT contains the in situ location!')
    
    return ofname