
class ComparisonOptions:

    def __init__(self):
        pass


    def hypstar_l2_spectral_variables(self):
        spectral_variables = ['water_leaving_radiance','reflectance','reflectance_nosc','std_water_leaving_radiance','std_reflectance','std_reflectance_nosc']
        return spectral_variables

    def hypstar_l1_spectral_variables(self):
        spectral_variables = ['downwelling_radiance','upwelling_radiance','irradiance','std_irradiance','std_downwelling_radiance','std_upwelling_radiance']
        return spectral_variables

    def hypstar_l2_single_variables(self):
        single_variables = ['rhof','n_total_scans','n_valid_scans','quality_flag','pointing_azimuth_angle','solar_azimuth_angle','solar_zenith_angle','viewing_azimuth_angle','viewing_zenith_angle']
        return single_variables

    def hypstar_l1_single_variables(self):
        single_variables = ['rhof_raa','rhof_sza','rhof_vza','rhof_wind']
        return single_variables

    def hypstar_aeronet_spectral_variables(self):
        import numpy as np

        spectral_variables = {
            'water_leaving_radiance':{
                'name':'Lw',
                'long_name': 'Water - leaving radiance',
                'units': 'mW/(cm2·μm·sr)',
                'scale_factor': 0.1
            },
            'downwelling_radiance':{
                'name': 'Li_mean',
                'long_name': 'Sky radiance - Mean',
                'units': 'mW/(cm2·μm·sr)',
                'scale_factor': 0.1
            },
            'upwelling_radiance':{
                'name': 'Lt_mean',
                'long_name': 'Total radiance from the sea surface - Mean',
                'units': 'mW/(cm2·μm·sr)',
                'scale_factor': 0.1

            },
            'irradiance':{
                'name': 'F0',
                'long_name': 'Extraterrestrial Solar Irradiance',
                'units': 'W/(m2·nm)',
                'scale_factor': 0.001
            },
            'reflectance':{
                'name': 'Rrs',
                'long_name': 'Reflectance of the water column at the surface',
                'units': 'sr-1',
                'scale_factor': 1/np.pi
            },
            'reflectance_nosc': {
                'name': 'Rrs_nosc',
                'long_name': 'Reflectance of the water column at the surface  without correction for the NIR similarity spectrum (see Ruddick et al., 2006)',
                'units': 'sr-1',
                'scale_factor': 1/np.pi
            },
            'std_downwelling_radiance':{
                'name': 'Li_stddev',
                'long_name': 'Sky radiance - Standard deviation',
                'units': 'mW/(cm2·μm·sr)',
                'scale_factor': 0.1
            },
            'std_upwelling_radiance': {
                'name': 'Lt_stddev',
                'long_name': 'Total radiance from the sea surface - Standard deviation',
                'units': 'mW/(cm2·μm·sr)',
                'scale_factor': 0.1
            }

        }

        return spectral_variables


    def aeronet_spectral_variables(self):
        spectral_variables = {
            'Rho':{
                'long_name': 'Sea-surface reflectance factor'
            },
            'Exact_Wavelengths': {
                'long_name': 'Exact wavelenghts',
                'units': 'nm'
            },
            'Lw':{
                'long_name': 'Water-leaving radiance',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Aerosol_Optical_Depth':{
                'long_name': 'Aerosol optical depth',
            },
            'F0': {
                'long_name': 'Extraterrestrial Solar Irradiance (F0) from Thuillier',
                'units': 'W/(m2·nm)'
            },
            'Li_mean': {
                'long_name': 'Sky radiance - Mean',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Li_stddev': {
                'long_name': 'Sky radiance - Standard deviation',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lt_mean': {
                'long_name': 'Total radiance from the sea surface - Mean',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lt_min_rel': {
                'long_name': 'Total radiance from the sea surface - Relative minimum',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lt_stddev': {
                'long_name': 'Total radiance from the sea surface - Standard deviation',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lwn': {
                'long_name': 'Normalized water-leaving radiance',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lwn_f_Q': {
                'long_name': 'Water-leaving radiance f/Q ',
                'units': 'mW/(cm2·μm·sr)'
            },
            'Lwn_IOP': {
                'long_name': 'Water-leaving radiance IOP',
                'units': 'mW/(cm2·μm·sr)'
            },
            'LwQ': {
                'long_name': 'Water-leaving radiance Q',
                'units': 'mW/(cm2·μm·sr)'
            },
            'NO2_Optical_Depth': {
                'long_name': 'NO2 Optical Depth',
            },
            'Ozone_Optical_Depth': {
                'long_name': 'Ozone Optical Depth',
            },
            'Rayleigh_Optical_Depth': {
                'long_name': 'Rayleigh Optical Depth',
            },
            'Water_Vapor_Optical_Depth': {
                'long_name': 'Water Varpor Optical Depth',
            },
            'Solar_Azimuth_Angle': {
                'long_name': 'Solar Azimuth Angle',
                'units': 'degrees'
            },
            'Solar_Zenith_Angle': {
                'long_name': 'Solar Zenith Angle',
                'units': 'degrees'
            }
        }

        return spectral_variables

    def aeronet_single_variables(self):
        single_variables = {
            'Chlorophyll-a': {
                'long_name': 'Chlorophyll-a concentration',
                'units': 'mg m-3'
            },
            'Number_of_Wavelengths':{
                'long_name': 'Number of wavelengths'
            },
            'Pressure(hPa)':{
                'long_name': 'Atmospheric pressure',
                'units': 'hPa'
            },
            'Total_NO2(DU)': {
                'long_name': 'Total NO2',
                'units': 'DU'
            },
            'Total_Ozone(Du)': {
                'long_name': 'Total Ozone',
                'units': 'Du'
            },
            'Total_Precipitable_Water(cm)': {
                'long_name': 'Total Precipitable Water',
                'units': 'cm'
            },
            'Wind_Speed(m_s)': {
                'long_name': 'Wind speed',
                'units': 'm/s'
            },
            'Observing_Azimuth_Angle':{
                'long_name': 'Observing azimuth angle',
                'units': 'degrees'
            },
            'Observing_Zenith_Angle': {
                'long_name': 'Observing zenith angle',
                'units': 'degrees'
            }
        }
        return single_variables

    def basic_mu_variables(self):

        basic_mu_variables = {
            'mu_wavelength': {
                'long_name': 'Match-up Wavelength',
                'units': 'nm'
            },
            'mu_day_id': {
                'long_name': 'Day Index'
            },
            'mu_AERONET_sequence_id': {
                'long_name': 'AERONET-OC sequence ID'
            },
            'mu_HYPSTAR_sequence_id': {
                'long_name': 'HYPSTAR sequence ID'
            },
            'mu_time_diff': {
                'long_name': 'Time difference between HYPSTAR and AERONET-OC acquisitions',
                'unit': 'seconds'
            }

        }

        return basic_mu_variables

    def mu_variables(self):
        mu_variables = ['Lw','F0','Li_mean','Li_stddev','Lt_mean','Lt_stddev','Rrs','Rrs_nosc']
        return mu_variables

    def mu_variables_keys(self):
        mu_variables_keys = ['AERONET','HYPSTAR_TO_AERONET']
        return mu_variables_keys