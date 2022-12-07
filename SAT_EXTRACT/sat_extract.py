import netCDF4
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
import numpy.ma as ma


class SatExtract:
    def __init__(self, ofname):

        self.FILE_CREATED = True
        try:
            self.EXTRACT = Dataset(ofname, 'w', format='NETCDF4')
        except PermissionError:
            print('Permission denied: ', ofname)
            self.FILE_CREATED = False

        self.geometry_variables = {
            'SZA': {
                'long_name': 'Sun Zenith Angle',
                'units': 'degrees'
            },
            'SAA': {
                'long_name': 'Sun Azimuth Angle',
                'units': 'degrees'
            },
            'OZA': {
                'long_name': 'Observation Zenith Angle',
                'units': 'degrees'
            },
            'OAA': {
                'long_name': 'Observation Azimuth Angle',
                'units': 'degrees'
            }
        }

    def set_global_attributes(self, at):
        # Atributes
        self.EXTRACT.creation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
        satellite = at['satellite']
        platform = at['platform']
        sensor = at['sensor']
        res_str = at['res']
        self.EXTRACT.satellite = satellite
        self.EXTRACT.platform = platform
        self.EXTRACT.sensor = sensor
        self.EXTRACT.description = f'{satellite}{platform} {sensor.upper()} {res_str} L2 extract'
        self.EXTRACT.satellite_aco_processor = at['aco_processor']  # 'Atmospheric Correction processor: xxx'
        self.EXTRACT.satellite_proc_version = at['proc_version']  # proc_version_str

        self.EXTRACT.insitu_site_name = at['station_name']
        self.EXTRACT.insitu_lat = at['in_situ_lat']
        self.EXTRACT.insitu_lon = at['in_situ_lon']

    def create_dimensions_basic(self, size_box):
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)

    def create_dimensions(self, size_box, n_bands):
        # dimensions
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)
        self.EXTRACT.createDimension('satellite_bands', n_bands)

    def create_dimensions_incluidinginsitu(self, size_box, n_bands, n_insitubands, n_insituid):
        # dimensions
        self.EXTRACT.createDimension('satellite_id', None)
        self.EXTRACT.createDimension('rows', size_box)
        self.EXTRACT.createDimension('columns', size_box)
        self.EXTRACT.createDimension('satellite_bands', n_bands)

        self.EXTRACT.createDimension('insitu_original_bands', n_insitubands)
        self.EXTRACT.createDimension('insitu_id', n_insituid)

    def create_geometry_variable(self, name_geom):
        if not name_geom in self.geometry_variables.keys():
            return None
        gvar = self.EXTRACT.createVariable(name_geom, 'f4', ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                           zlib=True, complevel=6)
        gvar.units = self.geometry_variables[name_geom]['units']
        gvar.long_name = self.geometry_variables[name_geom]['long_name']
        return gvar

    def create_satellite_time_variable(self, satellite_start_time):

        satellite_time = self.EXTRACT.createVariable('satellite_time', 'f8', ('satellite_id'), fill_value=-999,
                                                     zlib=True, complevel=6)
        #print('Satellite start time es: ',satellite_start_time)
        satellite_time[0] = float(satellite_start_time.timestamp())
        satellite_time.units = "Seconds since 1970-1-1"

    def create_pdu_variable(self, pdu, sensor):
        satellite_PDU = self.EXTRACT.createVariable('satellite_PDU', 'S4', ('satellite_id'), zlib=True,
                                                    complevel=6)  # string

        satellite_PDU[0] = pdu
        satellite_PDU.long_name = f'{sensor} source PDU name'

    def create_lat_long_variables(self, lat, lon, window):

        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        nrows = (window[1] - window[0])
        ncols = (window[3] - window[2])


        # latitude
        satellite_latitude = self.EXTRACT.createVariable('satellite_latitude', 'f8',
                                                         ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                         zlib=True, complevel=6)

        if lat.ndim == 1:
            for r in range(nrows):
                #satellite_latitude[0, r, :] = [lat[start_idx_x:stop_idx_x]]
                satellite_latitude[0, r, :] = [lat[start_idx_y:stop_idx_y]]
        else:
            satellite_latitude[0, :, :] = [lat[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
        satellite_latitude.short_name = 'latitude'

        # longitude
        satellite_longitude = self.EXTRACT.createVariable('satellite_longitude', 'f8',
                                                          ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                          zlib=True, complevel=6)
        if lon.ndim == 1:
            for c in range(ncols):
                #satellite_longitude[0, :, c] = [lon[start_idx_y:stop_idx_y]]
                satellite_longitude[0, :, c] = [lon[start_idx_x:stop_idx_x]]
        else:
            satellite_longitude[0, :, :] = [lon[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x]]
        satellite_longitude.short_name = 'longitude'

    def create_satellite_bands_variable(self, wavelengths):
        # Variable satellite_bands (wavelenghts for Rrs)
        satellite_bands = self.EXTRACT.createVariable('satellite_bands', 'f4', ('satellite_bands'), fill_value=-999,
                                                      zlib=True, complevel=6)
        satellite_bands[:] = wavelengths
        satellite_bands.units = 'nm'

    def create_rrs_variable(self, sensor):
        # Variable satellite_Rrs (NOT BRDF-corrected remote sensing reflectance)
        satellite_Rrs = self.EXTRACT.createVariable('satellite_Rrs', 'f4',
                                                    ('satellite_id', 'satellite_bands', 'rows', 'columns'),
                                                    fill_value=-999, zlib=True, complevel=6)
        satellite_Rrs.short_name = 'Satellite Rrs'
        satellite_Rrs.long_name = f'Above water Remote Sensing Reflectance for {sensor} acquisition'
        satellite_Rrs.units = "sr-1"
        return satellite_Rrs

    def create_aot_variable(self, aot, window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[3]
        stop_idx_x = window[4]
        # AOT
        satellite_AOT_0865p50_box = self.EXTRACT.createVariable('satellite_AOT_0865p50', 'f4',
                                                                ('satellite_id', 'rows', 'columns'), fill_value=-999,
                                                                zlib=True, complevel=6)
        satellite_AOT_0865p50_box[0, :, :] = ma.array(aot[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])
        satellite_AOT_0865p50_box.description = 'Satellite Aerosol optical thickness'

    def create_flag_variable(self, var_name, var_array, description, flag_masks, flag_meanings, window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        # Quality Flags
        satellite_flag = self.EXTRACT.createVariable(var_name, 'f4', ('satellite_id', 'rows', 'columns'),
                                                     fill_value=-999, zlib=True, complevel=6)
        satellite_flag[0, :, :] = [ma.array(var_array[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
        satellite_flag.description = description
        satellite_flag.flag_masks = flag_masks
        satellite_flag.flag_meanings = flag_meanings

    def create_2D_variable_general(self,var_name,var_array,window):
        start_idx_y = window[0]
        stop_idx_y = window[1]
        start_idx_x = window[2]
        stop_idx_x = window[3]
        dtype = var_array.datatype.str
        if dtype.startswith('<'):
            dtype = dtype[1:]
        satellite_2d_band = self.EXTRACT.createVariable(var_name,dtype,('satellite_id', 'rows', 'columns'),fill_value=-999, zlib=True, complevel=6)
        satellite_2d_band[0, :, :] = [ma.array(var_array[start_idx_y:stop_idx_y, start_idx_x:stop_idx_x])]
        return satellite_2d_band

    def create_insitu_time_variable(self):
        insitu_time = self.EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id', 'insitu_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'
        return insitu_time

    def create_insitu_original_bands_variable(self):
        insitu_original_bands = self.EXTRACT.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                            fill_value=-999, zlib=True, complevel=6)
        insitu_original_bands.description = 'In situ bands in nm'
        return insitu_original_bands

    def create_insitu_exact_wavelengths_variable(self):
        insitu_exact_wavelenghts = self.EXTRACT.createVariable('insitu_exact_wavelenghts', 'f4',
                                                               ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                               fill_value=-999, zlib=True, complevel=6)
        insitu_exact_wavelenghts.description = 'In situ bands in nm'

        return insitu_exact_wavelenghts

    def create_insitu_rrs_variable(self):
        insitu_Rrs = self.EXTRACT.createVariable('insitu_Rrs', 'f4',
                                                 ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                                 fill_value=-999, zlib=True, complevel=6)
        insitu_Rrs.description = 'In situ Rrs'

        return insitu_Rrs

    def create_insitu_time_difference_variable(self):

        time_difference = self.EXTRACT.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'),
                                                      fill_value=-999,
                                                      zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        return time_difference

    def create_insitu_lat_long_variables(self):

        insitu_lat = self.EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lat.short_name = "latitude"
        insitu_lat.units = "degrees"

        insitu_lon = self.EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lon.short_name = "longitude"
        insitu_lon.units = "degrees"



        return insitu_lat, insitu_lon

    def create_insitu_variables_for_single_insitu_data(self):

        insitu_lat = self.EXTRACT.createVariable('insitu_latitude', 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lat.short_name = "latitude"
        insitu_lat.units = "degrees"

        insitu_lon = self.EXTRACT.createVariable('insitu_longitude', 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_lon.short_name = "longitude"
        insitu_lon.units = "degrees"

        time_difference = self.EXTRACT.createVariable('time_difference', 'f4', ('satellite_id',),
                                                      fill_value=-999,
                                                      zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        insitu_time = self.EXTRACT.createVariable('insitu_time', 'f8', ('satellite_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        return insitu_lat,insitu_lon,insitu_time,time_difference

    def create_insitu_variable_for_single_insitu_data(self,name_var,units,desc):

        name = f'insitu_{name_var}'

        insitu_var = self.EXTRACT.createVariable(name, 'f8', ('satellite_id',),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        insitu_var.units = units
        insitu_var.short_name = name_var
        insitu_var.description = desc

        return insitu_var

    def create_insitu_flag_variable(self,name_var,flag_values,flag_meanings):

        name = f'insitu_flag_{name_var}'
        insitu_var = self.EXTRACT.createVariable(name,'i4', ('satellite_id',),fill_value=-999,zlib=True,complevel=6)
        insitu_var.short_name = name_var
        insitu_var.flag_masks = flag_values
        insitu_var.flag_meanings = flag_meanings

        return insitu_var

    def close_file(self):
        self.EXTRACT.close()
