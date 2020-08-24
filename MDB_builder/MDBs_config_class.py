class MDBs_config(object):
    """
    Configuration options for MDBs building. Attributes shall be edited before running MDBs builder.
    """
    def __init__(self):
        pass

    def defaults_olci(self):
        #comon paramters for Level 2
        import numpy as np
        import os
        from string import Template
        self.sensor = 'olci'
        self.satellite = 'S3'
        self.n_sensor_bands = 16 
        self.band_index = np.array([1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,21])
        self.file_size = 25
        self.product_template = {1: {'a':'S3A_OL_1_E','b':'S3B_OL_1_E'}, 2: {'a':'S3A_OL_2_W','b':'S3B_OL_2_W'} }
        self.var_ref = 'Oa01_reflectance' #used to check 
        #full path to ncf containing costant info for 16 bands
        self.data_dir = (os.environ['CODE_HOME']+'/DATA/')
        self.adf_OCP = {'A': self.data_dir + '/S3A_OL_2_OCP_AX_20160216T000000_20991231T235959_20170915T120000___________________MPC_O_AL_004.SEN3/OL_2_OCP_AX.nc','B':self.data_dir + '/S3B_OL_2_OCP_AX_20180425T000000_20991231T235959_20191023T120000___________________MPC_O_AL_002.SEN3/OL_2_OCP_AX.nc'}
        self.SRF = {'A':self.data_dir+'/S3A_OL_SRF_20160713_mean_rsr.nc4','B': self.data_dir+'/S3B_OL_SRF_0_20180109_mean_rsr.nc4'}
        self.sensor_bands = np.array ([400.0,412.5,442.5,490.0,510.0,560.0,620.0,665.0,673.75,681.25,708.75,753.75,778.75,865.0,885.0,1020.0])
        self.BRDF_bands = np.array([400.0,412.5,442.5,490.0,510.0,560.0,620.0,665.0,673.75,681.25,708.75])
        self.band_units = 'nm'
        self.reflectance_bands = [('Oa'+"{:02d}"+'_reflectance').format(curr) for curr in ([1,2,3,4,5,6,7,8,9,10,11,12,16,17,18,21])]
        self.sat_prefix = 'satellite_'
        self.insitu_prefix = 'insitu_'

        #this list is what's expected in miniprods to be complete and correct
        self.complete_list = (['time_stamp','Oa01_reflectance','Oa01_reflectance_err','Oa02_reflectance','Oa02_reflectance_err',
            'Oa03_reflectance','Oa03_reflectance_err','Oa04_reflectance','Oa04_reflectance_err','Oa05_reflectance',
            'Oa05_reflectance_err','Oa06_reflectance','Oa06_reflectance_err','Oa07_reflectance','Oa07_reflectance_err',
            'Oa08_reflectance','Oa08_reflectance_err','Oa09_reflectance','Oa09_reflectance_err','Oa10_reflectance',
            'Oa10_reflectance_err','Oa11_reflectance','Oa11_reflectance_err','Oa12_reflectance','Oa12_reflectance_err',
            'Oa16_reflectance','Oa16_reflectance_err','Oa17_reflectance','Oa17_reflectance_err','Oa18_reflectance',
            'Oa18_reflectance_err','Oa21_reflectance','Oa21_reflectance_err','CHL_OC4ME','CHL_OC4ME_err','CHL_NN','CHL_NN_err',
            'TSM_NN','TSM_NN_err','ADG443_NN','ADG443_NN_err','KD490_M07','KD490_M07_err','PAR','PAR_err','A865',
            'A865_err','T865','T865_err','WQSF','atmospheric_temperature_profile','horizontal_wind','humidity',
            'reference_pressure_level','sea_level_pressure','total_columnar_water_vapour','total_ozone','OAA','OZA','SAA','SZA',
            'FWHM','frame_offset','detector_index','lambda0','relative_spectral_covariance','solar_flux','latitude',
            'longitude','altitude','BRDF_Bands','BRDF'])
            
        ##NB: mminifies shall be organised in subdirectories as dir_structure/year/doy/
        self.miniprod_pattern = {'a': "S3A_OL_2*_25x25.nc", 'b': "S3B_OL_2*_25x25.nc"}

        #subdirectories structure used to look for miniprods
        self.dir_structure = {'a' : Template('/s3a/olci/wfr/nt/$label/'), 'b': Template('/s3b/olci/wfr/nt/$label/')}
        self.redo_extraction = True #True to redo corrupted or incmoplete extractions
        self.sub_dir_year_doy = True #set to false in case all the are in the same directory. Set to true to optimize research.
        self.file_exclude = 'NR'#timeliness to be escluded
        #time difference for matchup in MDBs
        self.time_diff = 86400 #seconds

        
        #OLCI variables to be included in MDBs, sorted as desired
        
        #reflectance are converted in Rrs corrected for BRDF before being ingested. 
        #time_stamp, Rrs and BRDF are required for MDB_reader
        self.first_variables = (['time_stamp','Oa01_reflectance',
                                 'Oa01_reflectance_err','Oa02_reflectance','Oa02_reflectance_err','Oa03_reflectance','Oa03_reflectance_err','Oa04_reflectance','Oa04_reflectance_err',
                                 'Oa05_reflectance','Oa05_reflectance_err','Oa06_reflectance','Oa06_reflectance_err','Oa07_reflectance','Oa07_reflectance_err','Oa08_reflectance','Oa08_reflectance_err',
                                 'Oa09_reflectance','Oa09_reflectance_err','Oa10_reflectance','Oa10_reflectance_err','Oa11_reflectance','Oa11_reflectance_err','Oa12_reflectance','Oa12_reflectance_err',
                                 'Oa16_reflectance','Oa16_reflectance_err','Oa17_reflectance','Oa17_reflectance_err', 'Oa18_reflectance','Oa18_reflectance_err','Oa21_reflectance','Oa21_reflectance_err','BRDF',
                                 'BRDF_Bands']) 
                
        #variables are divided and listed to force this ordering in the MDB        
        self.second_variables = (['CHL_OC4ME','CHL_OC4ME_err','CHL_NN','CHL_NN_err','TSM_NN','TSM_NN_err','ADG443_NN','ADG443_NN_err','KD490_M07','KD490_M07_err','PAR','PAR_err','A865',
                                  'A865_err','T865','T865_err','WQSF','FWHM','lambda0','frame_offset','horizontal_wind','humidity','reference_pressure_level','relative_spectral_covariance',
                                  'sea_level_pressure','total_columnar_water_vapour','total_ozone','OAA','OZA','SAA','SZA','detector_index','solar_flux','latitude','longitude','altitude'])
        
        self.exclude_fields = (['relative_spectral_covariance','reference_pressure_level','frame_offset','Oa01_reflectance_err','Oa02_reflectance_err','Oa03_reflectance_err',
                                'Oa04_reflectance_err','Oa05_reflectance_err','Oa06_reflectance_err','Oa07_reflectance_err','Oa08_reflectance_err','Oa09_reflectance_err','Oa10_reflectance_err','Oa11_reflectance_err',
                                'Oa12_reflectance_err','Oa16_reflectance_err','Oa17_reflectance_err','Oa18_reflectance_err','Oa21_reflectance_err','CHL_OC4ME_err','CHL_NN_err','TSM_NN_err','ADG443_NN_err',
                                'KD490_M07_err','T865_err','A865_err','PAR_err'])
        self.append_fields = list() 
        
    def defaults_AERONET(self): 
        '''
        Specifies options for AERONET-OC data in MDBs. If any new AERONET site becomes available, it needs to be added, together with refeerences (PI and contacts)
        '''
        import numpy as np
        #all possible wavelengths in original data
        self.all_original_bands = np.array([411,412,413,440,441,442,488,489,490,491,530,550,551,552,553,554,555,666,667,668,868,869,870,871,1017,1018,1019,1020,1021])   
        #AERONET sites on lakes
        self.lakes = (['Palgrunden','Lake_Okeechobee','South_Greenbay','Lake_Erie'])
        #list of variable included in AERONET QA ncf files 
        self.aero_original_1d = (['Julian_Day','Instrument_Number','Pressure','Wind_Speed','Chlorophyll','Sea_Surface_Reflectance','Ozone','time','QA','SAM'])#('level','lat',
                                #'long','altitude' are attribute of the in situ ncf file and appended as costant values by
                                #pass_aeronet_dataset.py)
        #central nominal wavelength at which AERONET Version 2 data are grouped or shifted for Rrs. 
        self.Rrs_aeronet_bands = np.array([412.5,442.5,490.0,510.0,530.0,550.0,555.0,560.0,665.0,668.0,865.0,1020.0])          
        self.aero_Nbands = self.Rrs_aeronet_bands.size
        
        #list of AERONET variables to be included in MDBs, divided according to final dimensions and sorted as desired,
        #corresponding units and descriptions. Remeber that Lwn_fQ is stored in group in AERONET files
        self.groups3D = np.array(['Lw','Lw_Q','Lwn','Lwn_fQ','Lt_Mean','Lt_Stddev','Lt_Min_rel','Li_Mean_val','Li_Stddev',
                                 'AOT','OOT','ROT','Solar_Zenith','Solar_Azimuth','ExactWavelength','F0'])
        self.units3D = np.array(['mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1',
                                'mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1',
                                'mW.cm-2.sr-1.um-1','-','-','-','degrees','degrees','nm','mW.m-2.nm-1'])
        self.description3D = np.array(['AERONET-OC Water-Leaving Radiance','AERONET-OC Water-Leaving Radiance corrected for viewing angle dependence' ,'AERONET-OC Normalized Water-Leaving Radiance determined from Lw (i.e., not corrected for the viewing angle dependence and for the effects of the non-isotropic distribution of the in-water light field','AERONET-OC Normalized Water-Leaving Radiance determined from Lw_Q corrected for the effects of the non-isotropic distribution of the in-water radiance field (i.e., f/Q Corrected)', 'AERONET-OC Mean of 11 Above-Water Total Radiance Measurements','AERONET-OC Standard Deviation of 11 Above-Water Total Radiance Measurements','AERONET-OC Mean of Lowest Two of 11 Above-Water Total Radiance Measurements','AERONET-OC Mean of 3 Sky Radiance Measurements','AERONET-OC Standard Deviation of 3 Sky Radiance Measurement','AERONET-OC Aerosol Optical Thickness','AERONET-OC Ozone Optical Thickness','AERONET-OC Rayleigh Optical Thickness','AERONET-OC Solar Zenith Angle','AERONET-OC Solar Azimuth Angle','AERONET-OC exact wavelength','Thuillier et al. (2003) solar irradiance resampled for AERONET-OC exact wavelength, used for Lwn_fQ to Rrs conversion'])
        
        self.var2D_1 = np.array(['Julian_Day'])
        self.units2D_1 = np.array(['-'])
        self.description2D_1 = np.array(['AERONET-OC acquisition time in days from January 1st'])
        for b in self.band_index:
            self.var2D_1 = np.append(self.var2D_1, "Oa%02d_Rrs" %b) #mandatory
            self.description2D_1 = np.append(self.description2D_1,'Above water Remote Sensing Reflectance for AERONET-OC acquisition at %d nm obtained from Lwn_fQ' %float(self.sensor_bands[self.band_index == b]))
            self.units2D_1 = np.append(self.units2D_1,'sr-1')
        for b in self.band_index:
            self.var2D_1 = np.append(self.var2D_1,'Oa%02d_Rrs_applied_shift' %b) #mandatory 
            self.description2D_1 = np.append(self.description2D_1,'True if values at %d nm is derived trhough band_shifting correction' %float(self.sensor_bands[self.band_index == b]))
            self.units2D_1 = np.append(self.units2D_1,'-')
        self.units2D_1 = np.append(self.units2D_1,np.array(['mg.m-3','-','-','degrees']))
        self.description2D_1 = np.append(
                self.description2D_1,np.array(['AERONET-OC Chlorophyll-a','AERONET-OC Sea Surface Reflectance',
                                               'Quality level based on SAM','Spectral Angle Mapper for AERONET-OC spectrum']))
        self.var2D_1 = np.append(self.var2D_1,np.array(
            ['Chlorophyll','Sea_Surface_Reflectance','QA','SAM']))
        self.var2D_2 = np.array(['Instrument_Number','level','lat','long','altitude','Pressure','Wind_Speed','Ozone'])
        self.units2D_2 = np.array(['-','-','degrees_north','degrees_east','m','hPa','m.s-1','Dobson_Units'])
        self.description2D_2 = np.array(['AERONET-OC Instrument Number','AERONET-OC quality level',
                                         'Latitude for AERONET-OC site','Longitude for AERONET-OC site',
                                         'Altitude for AERONET-OC site','Atmospheric Pressure at AERONET-OC site',
                                         'Wind Speed at AERONET-OC site','Ozone at AERONET-OC site'])
        self.site_list = (['ARIAKE_TOWER','Blyth_NOAH','Casablanca_Platform','COVE_SEAPRISM','Galata_Platform','Gloria','GOT_Seaprism','Grizzly_Bay','Gustav_Dalen_Tower','HBOI','Helsinki_Lighthouse','Ieodo_Station','Irbe_Lighthouse','Kemigawa_Offshore',
        'Lake_Erie','Lake_Okeechobee','LISCO','Lucinda','MVCO','Palgrunden','Section-7_Platform','Socheongcho','South_Greenbay','Thornton_C-power','USC_SEAPRISM','USC_SEAPRISM_2','Venise','WaveCIS_Site_CSI_6','Zeebrugge-MOW1'])        
        self.PIs = {'ARIAKE_TOWER':'Joji Ishizaka, Kohei Arai',
                    'Bahia_Blanca':'Brent Holben',
                    'Blyth_NOAH':'Rodney Forster',
                    'Casablanca_Platform':'Giuseppe Zibordi, Marco Talone',
                    'COVE_SEAPRISM':'Brent Holben',
                    'Galata_Platform':'Giuseppe Zibordi',
                    'Gloria':'Giuseppe Zibordi',
                    'GOT_Seaprism':'Brent Holben',
                    'Grizzly_Bay':'Nima Pahlevan',
                    'Gustav_Dalen_Tower': 'Giuseppe Zibordi',
                    'HBOI':'Nima Pahlevan',
                    'Helsinki_Lighthouse':'Giuseppe Zibordi',
                    'Ieodo_Station':'Young-Je Park, Hak-Yeol You', 
                    'Irbe_Lighthouse':'Giuseppe Zibordi',
                    'Kemigawa_Offshore': 'Brent_Holben',
                    'Lake_Erie':'Tim Moore, Steve Ruberg, Menghua Wang',
                    'Lake_Okeechobee':'Nima Pahlevan',
                    'LISCO':'Sam Ahmed, Alex Gilerson',
                    'Lucinda':'Thomas Schroeder',
                    'MVCO':'Hui Feng, Heidi M. Sosik',
                    'Palgrunden':'Susanne Kratzer',
                    'Section-7_Platform':'Giuseppe Zibordi',
                    'Socheongcho':'Young-Je Park',
                    'South_Greenbay':'Nima Pahlevan',
                    'Thornton_C-power':'Dimitry Van der Zande',
                    'USC_SEAPRISM':'Nick Tufillaro',
                    'USC_SEAPRISM_2':'Nick Tufillaro',
                    'Venise':'Giuseppe Zibordi',
                    'WaveCIS_Site_CSI_6':'Alan Weidemann, Bill Gibson, Robert Arnone',
                    'Zeebrugge-MOW1':'Dimitry Van der Zande'                  
                    }
        self.PIs_mail = {'ARIAKE_TOWER':'Jishizaka@nagoya-u.jp, arai@cc.saga-u.ac.jp',
                         'Blyth_NOAH':'r.forster@hull.ac.uk',
                         'Bahia_Blanca':'Brent.N.Holben@nasa.gov',
                         'Casablanca_Platform':'giuseppe.zibordi@ec.europa.eu, marco.talone@ec.europa.eu',
                         'COVE_SEAPRISM':'Brent.N.Holben@nasa.gov',
                         'Galata_Platform':'giuseppe.zibordi@ec.europa.eu',
                         'Gloria':'giuseppe.zibordi@ec.europa.eu',                         
                         'GOT_Seaprism':'Brent.N.Holben@nasa.gov',
                         'Grizzly_Bay': 'nima.pahlevan@nasa.gov',
                         'Gustav_Dalen_Tower':'giuseppe.zibordi@ec.europa.eu',
                         'HBOI':'nima.pahlevan@nasa.gov',
                         'Helsinki_Lighthouse':'giuseppe.zibordi@ec.europa.eu',
                         'Ieodo_Station':'youngjepark@kiost.ac,peterhak@korea.kr',
                         'Irbe_Lighthouse':'giuseppe.zibordi@ec.europa.eu',
                         'Kemigawa_Offshore':'Brent.N.Holben@nasa.gov',
                         'Lake_Erie':'timothy.moore@unh.edu, steve.ruberg@noaa.gov, Menghua.Wang@noaa.gov',
                         'Lake_Okeechobee':'nima.pahlevan@nasa.gov',
                         'LISCO':'ahmed@ccny.cuny.edu, gilerson@ccny.cuny.edu',
                         'Lucinda':'Thomas.Schroeder@csiro.au',
                         'MVCO':'Hui.Feng@unh.edu, hsosik@whoi.edu',
                         'Palgrunden':'Susanne.Kratzer@su.se',
                         'Section-7_Platform':'giuseppe.zibordi@ec.europa.eu',
                         'Socheongcho':'youngjepark@kiost.ac',
                         'South_Greenbay':'nima.pahlevan@nasa.gov',
                         'Thornton_C-power':'dvanderzande@naturalsciences.be',
                         'USC_SEAPRISM':'nbt@coas.oregonstate.edu',
                         'USC_SEAPRISM_2':'nbt@coas.oregonstate.edu',
                         'Venise':'giuseppe.zibordi@ec.europa.eu',
                         'WaveCIS_Site_CSI_6':'Alan.Weidemann@nrlssc.navy.mil, bgibson@lsu.edu, Robert.Arnone@usm.edu',
                         'Zeebrugge-MOW1':'dvanderzande@naturalsciences.be'
                         }

        #AERONET MDB attributes
        self.description = "S3A OLCI WFR L2 - AERONET-OC Matchups Data Base"
        self.description_B = "S3B OLCI WFR L2 - AERONET-OC Matchups Data Base"
        self.description_L1 = "S3A OLCI EFR L1 - AERONET-OC Matchups Data Base"
        self.description_L1_B = "S3B OLCI EFR L1 - AERONET-OC Matchups Data Base"
        self.data_usage = ('Notice to users: AERONET-OC data have been downloaded from https://aeronet.gsfc.nasa.gov/cgi-bin/type_piece_of_map_seaprism_new. The applicable data policies must be followed. Please check at https://aeronet.gsfc.nasa.gov/cgi-bin/type_piece_of_map_seaprism_new about Data Usage Policy')
        
        #list of variables to be excluded
        self.aero_exclude = (['QA'])
        
        ##params for banshifting
        self.QAA_g0 = 0.08945 #Melin&Sclep, 2015
        self.QAA_g1 = 0.1247   #Melin&Sclep, 2015 
        self.DlambdaMAX_vis = 1 #nm: max distance for which bandshifting is not required in VIS
        self.DlambdaMAX_NIR = 6 #nm: max distance for which bands are compared in NIR
        
    def defaults_aeronet_V3(self):
        '''
        Updates optios for AERONET-OC data in MDBs when using AERONET-OC version 3 data
        '''        
        import numpy as np
        self.Rrs_aeronet_bands = self.sensor_bands
        self.aero_Nbands = self.Rrs_aeronet_bands.size
        #list of variable included in AERONET QA ncf files version 3
        self.aero_original_1d = (['Day_of_Year_fraction','AERONET_Instrument_Number','Pressure','Wind_Speed','Chlorophyll','Total_Ozone','time','SAM'])#,'QA'])
        self.var2D_1 = np.array(['Day_of_Year_fraction'])
        self.units2D_1 = np.array(['-'])
        self.description2D_1 = np.array(['AERONET-OC acquisition time in days from January 1st'])
        for b in self.band_index:
            self.var2D_1 = np.append(self.var2D_1, "Oa%02d_Rrs" %b) #mandatory
            self.description2D_1 = np.append(self.description2D_1,'Above water Remote Sensing Reflectance for AERONET-OC acquisition at %d nm obtained from Lwn_fQ' %float(self.sensor_bands[self.band_index == b]))
            self.units2D_1 = np.append(self.units2D_1,'sr-1')
        for b in self.band_index:
            self.var2D_1 = np.append(self.var2D_1,'Oa%02d_Rrs_applied_shift' %b) #mandatory 
            self.description2D_1 = np.append(self.description2D_1,'True if values at %d nm is derived trhough band_shifting correction' %float(self.sensor_bands[self.band_index == b]))
            self.units2D_1 = np.append(self.units2D_1,'-')
        self.units2D_1 = np.append(self.units2D_1,np.array(['mg.m-3','degrees']))#'-',]))
        self.description2D_1 = np.append(self.description2D_1,np.array(['AERONET-OC Chlorophyll-a','Spectral Angle Mapper for AERONET-OC spectrum']))#,'Quality level based on SAM']))
        self.var2D_1 = np.append(self.var2D_1,np.array(['Chlorophyll','SAM'])) #,'QA']))
        self.var2D_2 = np.array(['AERONET_Instrument_Number','level','lat','long','altitude','Pressure','Wind_Speed','Total_Ozone'])
        self.units2D_2 = np.array(['-','-','degrees_north','degrees_east','m','hPa','m.s-1','Dobson_Units'])
        self.description2D_2 = np.array(['AERONET-OC Instrument Number','AERONET-OC quality level','Latitude for AERONET-OC site','Longitude for AERONET-OC site',
                                         'Altitude for AERONET-OC site','Atmospheric Pressure at AERONET-OC site','Wind Speed at AERONET-OC site','Ozone at AERONET-OC site'])
        self.groups3D = np.array(['Lw','Lw_Q','Lwn','Lwn_fQ','Lt_mean','Lt_stddev','Lt_min_rel','Li_mean','Li_stddev','Rho',
                                 'Aerosol_Optical_Depth','Ozone_Optical_Depth','Rayleigh_Optical_Depth','Solar_Zenith_Angle',
                                  'Solar_Azimuth_Angle','Exact_Wavelengths','F0'])
        self.units3D = np.array(['mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1',
                                'mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1','mW.cm-2.sr-1.um-1',
                                'mW.cm-2.sr-1.um-1','-','-','-','-','degrees','degrees','nm','mW.m-2.nm-1'])
        self.description3D = np.array(['AERONET-OC Water-Leaving Radiance','AERONET-OC Water-Leaving Radiance corrected for viewing angle dependence' ,'AERONET-OC Normalized Water-Leaving Radiance determined from Lw (i.e., not corrected for the viewing angle dependence and for the effects of the non-isotropic distribution of the in-water light field','AERONET-OC Normalized Water-Leaving Radiance determined from Lw_Q corrected for the effects of the non-isotropic distribution of the in-water radiance field (i.e., f/Q Corrected)', 'AERONET-OC Mean of 11 Above-Water Total Radiance Measurements','AERONET-OC Standard Deviation of 11 Above-Water Total Radiance Measurements','AERONET-OC Mean of Lowest Two of 11 Above-Water Total Radiance Measurements','AERONET-OC Mean of 3 Sky Radiance Measurements','AERONET-OC Standard Deviation of 3 Sky Radiance Measurement','AERONET-OC Sea Surface Reflectance','AERONET-OC Aerosol Optical Thickness','AERONET-OC Ozone Optical Thickness','AERONET-OC Rayleigh Optical Thickness','AERONET-OC Solar Zenith Angle','AERONET-OC Solar Azimuth Angle','AERONET-OC exact wavelength','Thuillier et al. (2003) solar irradiance resampled for AERONET-OC exact wavelength, used for Lwn_fQ to Rrs conversion'])
        self.data_usage = ('Notice to users: AERONET-OC data have been downloaded from https://aeronet.gsfc.nasa.gov/cgi-bin/draw_map_display_seaprism_v3. The applicable data policies must be followed. Please check at https://aeronet.gsfc.nasa.gov/cgi-bin/draw_map_display_seaprism_v3 about Data Usage Policy')
        
    def defaults_MOBY(self): 
        import numpy as np
        #variabless to be included in MDBs
        self.moby_variables = (['Rrs','Lwn','Es']) #Rrs required! (for MDB_reader). If Lwn included , Lw is included by default
        self.moby_Rrs_bands = [('Oa'+"{:02d}"+'_Rrs').format(curr) for curr in range(1,22)]
        self.moby_wave = ([400.0,412.5,442.5,490.0,510.0,560.0,620.0,665.0,673.75,681.25,708.75,753.75,761.25,764.375,767.5,778.75,865.0,885.0,900.0,940.0,1020.0])
        
        self.MOBY_no_24 = 60 #no. of possible measurements available within +-24 hours from overpass. Up to now is 12 but may increase with new sensors launches
        #moby Es bands to be included
        
        #subdiretory structure for MOBY files
        self.MOBY_insitu_dir_struct = {'a':'/MOBY/S3A/OLCI/', 'b': '/MOBY/S3A/OLCI/'}
        #deployments with post-deployment calibration applied (check on MOBY site)
        self.MOBY_list_postCal = (259,260,261,262,263,264,266,267) #post cals is available for deployment up to 267
        
        #additional variables to be included in MOBY MDB. 
        self.additional_variables = np.array(['Year','JDay','Deploy','Lat','Long','Data','GMTTime','arm','ObsDate','process_status'])
        self.additional_description = (['Year of the corresponding MOBY acquisition',
                                      'Day of the year of the corresponding MOBY acquisition','MOBY deployment identification number',
                                      'Latitude for MOBY buoy','Longitude for MOBY buoy',
                                      'MOBY data status (1 = Good and 2 = Questionable)','GMTTime','MOBY arms used: 1 = LuTop-LuMid, 2 = LuTop-LuBot, 7 = LuMid-LuBot', 'ObsDate', 'MOBY data reprocessing status: pre (0) or post (1) deployment calibration'])  
        self.additional_units = (['years','Days since November 24th, 4714 BC at 12:00 in the proleptic Gregorian calendar','-','degrees_north','degrees_east','-','-','-','-','-'])
        #list of variable to be excluded
        self.moby_exclude = ()
        #MDB attributes
        self.MOBYdescription = {'A':"S3A OLCI WFR L2 - MOBY Matchups Data Base", 'B': "S3B OLCI WFR L2 - MOBY Matchups Data Base"}
        self.MOBYdescription_L1 = {'A':"S3A OLCI EFR L1 - MOBY Matchups Data Base", 'B':"S3B OLCI EFR L1 - MOBY Matchups Data Base"}
        self.MOBYdata_usage = ('MOBY data have been downloaded from https://www.star.nesdis.noaa.gov/sod/moby/gold/. '
                               'We acknowledge NOAA CoastWatch/OceanWatch for them. The applicable data policies must be followed. Please check at https://www.mlml.calstate.edu/moby/ about data usage policy.') 
        

    def defaults_cruises(self): 
        '''
        Configuration options for in situ data in MDBs when using data different that AERONET-OC or MOBY
        ''' 
        #variables in sources files
        self.cruises_variables = (['time','identifier_product_doi','received','source','investigators','affiliations','contact','experiment','cruise','data_file_name','data_type','data_status','cloud_percent','wave_height','wind_speed','secchi_depth','water_depth','lat','long','HPLC','CHL','min_depth','ADG443','min_depth_ad','min_depth_ag','KD490','KD490_sd','Es_stability','Rrs_400.0','Rrs_412.5','Rrs_442.5','Rrs_490.0','Rrs_510.0','Rrs_560.0','Rrs_620.0','Rrs_665.0','Rrs_673.75','Rrs_681.25','Rrs_708.75','Rrs_753.75','Rrs_778.75','Rrs_865.0','Rrs_885.0','Rrs_1020.0','Rrs_400.0_sd','Rrs_412.5_sd','Rrs_442.5_sd','Rrs_490.0_sd','Rrs_510.0_sd','Rrs_560.0_sd','Rrs_620.0_sd','Rrs_665.0_sd','Rrs_673.75_sd','Rrs_681.25_sd','Rrs_708.75_sd','Rrs_753.75_sd','Rrs_778.75_sd','Rrs_865.0_sd','Rrs_885.0_sd','Rrs_1020.0_sd','QA'])
        self.cruises_description = (['In situ measurement acquisition time','Collection DOI','Date of data submission by PI(s)','Source Database name','Cruise investigator(s)','Cruise Investigator(s) affiliation','Cruise investigator(s) contact','Name of the experiment(s)','Cruise name','Source file(s) name', 'Data type','Data status','Cloud percentage observed in situ','Waves height observed in situ','Wind speed measured in situ','Secchi depth estimated in situ','Water depth measured in situ','In situ measurement latitude','In situ measurement longitude','Chlorophyll-a concentration from HPLC','Chlorophyll a derived fluorometrically/spectrophotometrically','Depth of the shallowest in situ measurement','In situ CDM absorption coefficient at 443 nm','Depth of the shallowest in situ measurement of ad','Depth of the shallowest in situ measurement of ag','In situ diffuse attenuation coefficient at 490 nm','In situ diffuse attenuation coefficient at 490 nm standard deviation','Flag for in situ irradiance stability','In situ Remote Sensing reflectance at 400.0 nm','In situ Remote Sensing reflectance at 412.5 nm','In situ Remote Sensing reflectance at 442.5 nm','In situ Remote Sensing reflectance at 490.0 nm','In situ Remote Sensing reflectance at 510.0 nm','In situ Remote Sensing reflectance at 560.0 nm','In situ Remote Sensing reflectance at 620.0 nm','In situ Remote Sensing reflectance at 665.0 nm','In situ Remote Sensing reflectance at 673.75 nm','In situ Remote Sensing reflectance at 682.15 nm','In situ Remote Sensing reflectance at 708.75 nm','In situ Remote Sensing reflectance at 753.75 nm','In situ Remote Sensing reflectance at 778.75 nm','In situ Remote Sensing reflectance at 865.0 nm','In situ Remote Sensing reflectance at 885.0 nm','In situ Remote Sensing reflectance at 1020.0 nm', 'In situ Remote Sensing reflectance at 400.0 nm standard deviation','In situ Remote Sensing reflectance at 412.5 nm standard deviation','In situ Remote Sensing reflectance at 442.5 nm standard deviation','In situ Remote Sensing reflectance at 490.0 nm standard deviation','In situ Remote Sensing reflectance at 510.0 nm standard deviation','In situ Remote Sensing reflectance at 560.0 nm standard deviation','In situ Remote Sensing reflectance at 620.0 nm standard deviation','In situ Remote Sensing reflectance at 665.0 nm standard deviation','In situ Remote Sensing reflectance at 673.75 nm standard deviation','In situ Remote Sensing reflectance at 682.15 nm standard deviation','In situ Remote Sensing reflectance at 708.75 nm standard deviation','In situ Remote Sensing reflectance at 753.75 nm standard deviation','In situ Remote Sensing reflectance at 778.75 nm standard deviation','In situ Remote Sensing reflectance at 865.0 nm standard deviation','In situ Remote Sensing reflectance at 885.0 nm standard deviation','In situ Remote Sensing reflectance at 1020.0 nm standard deviation','Quality flags for in situ Remote Sensing reflectance'])
        self.cruises_units = (['Seconds since 1970-1-1','-','-','-','-','-','-','-','-','-','-','-','-','m','m.s-1','m','m','degrees','degrees','mg.m-3','mg.m-3','m','m-1','m','m','m-1','m-1','-','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','sr-1','-'])
        #subdirectory structure for miniprods main directory
        #self.subdirectory_structure = '/S3A/OLCI/WFR/NATIVE/INST/EXTRACTON/IC01/2/'
        #self.subdirectory_structure_B = '/S3B/OLCI/WFR/NATIVE/INST/EXTRACTON/IC01/2/'
        
        #names of csv files containing Seabass data collected and QCed
        self.HPLC = 'HPLC.csv'
        self.CHL = 'CHL.csv'
        self.RRS = 'RRS.csv'
        self.KD490 = 'KD490.csv'
        self.ADG443 = 'ADG443.csv'
        self.TSM = 'TSM.csv'
        
        #maximum distance from scene central pixel and measuremnt station for matchup definition
        self.distance = 0.150 #km
        
        #thresholds for quality scores (by Lee et al., see INSITU_QA) to assign QA level to insitu RRS spectra
        self.highest_quality_RRS = 0.8
        self.less_quality_RRS = 0.6
        
        #list of variables to be excluded form MDBs
        self.insitu_exclude = (['date_EUM_proc', 'public'])
        
        #MDBs attributes
        self.insitu_description = "S3A OLCI WFR L2 - xx Matchups Data Base"
        self.insitu_description_B = "S3B OLCI WFR L2 - xx Matchups Data Base" # xx is replaced by var type
        self.insitu_description_L1 = "S3A OLCI EFR L1 - xx Matchups Data Base"
        self.insitu_description_L1_B = "S3B OLCI EFR L1 - xx Matchups Data Base"
        self.data_usage = ('Depending on the in situ data source, the applicable data policies shall be followed. Applicable data policy is available at https://ocdb.eumetsat.int/') 
        
    def defaults_common_L1(self):
        '''
        Updates optios  satellite products data in MDBs when using Level 1 products
        '''         
        import numpy as np
        from string import Template
        self.n_sensor_bands = 21
        self.file_size = 25
        self.var_ref = 'Oa01_radiance' #variable used to check minifile
        
        self.sensor_bands = np.array ([400.0,412.5,442.5,490.0,510.0,560.0,620.0,665.0,673.75,681.25,708.75,753.75,761.25,764.375,767.5,778.75,865.0,885.0,900.0,940.0,1020.0])
        self.band_index = np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21])
        #this list is what's expected in miniprods to be complete and correct
        self.complete_list = (['time_stamp','Oa01_radiance','Oa02_radiance','Oa03_radiance','Oa04_radiance','Oa05_radiance',
                               'Oa06_radiance','Oa07_radiance','Oa08_radiance','Oa09_radiance','Oa10_radiance','Oa11_radiance',
                               'Oa12_radiance','Oa13_radiance','Oa14_radiance','Oa15_radiance','Oa16_radiance','Oa17_radiance',
                               'Oa18_radiance','Oa19_radiance','Oa20_radiance','Oa21_radiance','quality_flags',
                               'atmospheric_temperature_profile','horizontal_wind','humidity','reference_pressure_level',
                               'sea_level_pressure','total_columnar_water_vapour','total_ozone','OAA','OZA','SAA','SZA',
                               'FWHM','frame_offset','detector_index','lambda0','relative_spectral_covariance','solar_flux',
                               'latitude','longitude','altitude'])

        #subdirectories structure used to look for miniprods
        self.miniprod_pattern = {'a': "S3A_OL_1*_25x25.nc", 'b' : "S3B_OL_1*_25x25.nc"}
        self.dir_structure = {'a' : Template('/s3a/olci/efr/nt/$label/'), 'b': Template('/s3b/olci/efr/nt/$label/')}

        #OLCI variables to be included in MDBs, sorted as desired
        #time_stamp, is required for MDB_reader
        self.first_variables =(['time_stamp','Oa01_radiance','Oa02_radiance','Oa03_radiance','Oa04_radiance','Oa05_radiance',
                               'Oa06_radiance','Oa07_radiance','Oa08_radiance','Oa09_radiance','Oa10_radiance','Oa11_radiance',
                               'Oa12_radiance','Oa13_radiance','Oa14_radiance','Oa15_radiance','Oa16_radiance','Oa17_radiance',
                               'Oa18_radiance','Oa19_radiance','Oa20_radiance','Oa21_radiance','quality_flags',
                               'atmospheric_temperature_profile','horizontal_wind','humidity','reference_pressure_level',
                               'sea_level_pressure','total_columnar_water_vapour','total_ozone','OAA','OZA','SAA','SZA',
                               'FWHM','frame_offset','detector_index','lambda0','relative_spectral_covariance','solar_flux',
                               'latitude','longitude','altitude'])
        self.second_variables = ([])
        self.exclude_fields = (['relative_spectral_covariance'])
        self.append_fields = list()
