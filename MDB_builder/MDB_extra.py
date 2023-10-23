import os
from netCDF4 import Dataset
from datetime import datetime as dt
import numpy.ma as ma


class MDBExtra():
    def __init__(self, path, verbose):
        self.info_path = self.extract_info_path(path, verbose)
        self.time_max = 120
        self.variables = None

    def extract_info_path(self, path, verbose):
        info_path = {}
        for name in os.listdir(path):
            if not name.endswith('.nc'):
                continue
            fpath = os.path.join(path, name)
            try:
                nc_sat = Dataset(fpath)
            except:
                continue
            if 'satellite' not in nc_sat.ncattrs() and 'platform' not in nc_sat.ncattrs() and 'satellite_time' not in nc_sat.variables:
                continue

            key = f'{nc_sat.satellite}{nc_sat.platform}'
            datehere = dt.fromtimestamp(float(nc_sat.variables['satellite_time'][0]))
            if verbose:
                print(f'[INFO] Name file: {name} Key: {key} Date: {datehere}')
            date_str = datehere.strftime('%Y%m%d')
            time_str = datehere.strftime('%H%M%S')
            if key not in info_path.keys():
                info_path[key] = {
                    date_str: {
                        'Time': time_str,
                        'Path': fpath
                    }
                }
            else:
                info_path[key][date_str] = {
                    'Time': time_str,
                    'Path': fpath
                }

        return info_path

    def get_extradaset(self, satellite, platform, sat_time, lat_array, lon_array):
        key = f'{satellite}{platform}'
        date_str = sat_time.strftime('%Y%m%d')
        #print(f'{satellite}{platform}{sat_time}========================================================================')
        ncdataset = None
        if key in self.info_path.keys():
            if date_str in self.info_path[key].keys():
                # print(self.info_path[key][date_str]['Time'])
                time_str = self.info_path[key][date_str]['Time']
                extradatestr = f'{date_str}T{time_str}'
                extradate = dt.strptime(extradatestr, '%Y%m%dT%H%M%S')
                valid = True
                if self.time_max > 0:
                    time_dif = abs((sat_time - extradate).total_seconds())
                    if time_dif > self.time_max:
                        print('TIME DIF > TIME MAX',key, sat_time, extradate, time_dif)
                        valid = False

                if valid:
                    ncdataset = Dataset(self.info_path[key][date_str]['Path'])
                    if lat_array is not None and lon_array is not None:

                        lat_array_here = ma.array(ncdataset.variables['satellite_latitude'][:])
                        lon_array_here = ma.array(ncdataset.variables['satellite_longitude'][:])

                        max_dif_lat = ma.max(ma.abs(lat_array_here - lat_array))
                        max_dif_lon = ma.max(ma.abs(lon_array_here - lon_array))
                        print(f'[INFO] Max dif lat: {max_dif_lat} Max dif lon: {max_dif_lon}')

                        #print(max_dif_lat,'<------->',max_dif_lon)

                        # if max_dif_lat > 0.000001 or max_dif_lon > 0.000001:
                        # #if max_dif_lat > 0.01 or max_dif_lon > 0.01:
                        #     ncdataset = None
                        #     print(f'[WARNING] Lat and long arrays of main and secondary extract do not coincide. Max dif lat: {max_dif_lat} Max dif lon: {max_dif_lon}')
                else:
                    print(
                        f'[WARNING] Extra sat extract are not available withing the corresponding time window for: {key} {sat_time}')
            else:
                print(f'[WARNING] Extra sat extracts are not available for date:  {key} {sat_time}')
        return ncdataset


def check(path_main, path_secondary):
    mdbe = MDBExtra(path_main, True)

    ntotal = 0
    nvalid = 0
    for name in os.listdir(path_secondary):
        if not name.endswith('.nc'):
            continue
        fname = os.path.join(path_secondary, name)
        ncsat = Dataset(fname)
        ntotal = ntotal + 1
        sat_time = dt.fromtimestamp(float(ncsat.variables['satellite_time'][0]))
        lat_array = ma.array(ncsat.variables['satellite_latitude'][:])
        lon_array = ma.array(ncsat.variables['satellite_longitude'][:])
        ncsatextra = mdbe.get_extradaset(ncsat.satellite, ncsat.platform, sat_time, lat_array, lon_array)


        if ncsatextra is not None:
            nvalid = nvalid + 1
            var_secondary = ncsatextra.variables['satellite_pixel_classif_flags']
            dtype = var_secondary.datatype.str
            if dtype.startswith('<'):
                dtype = dtype[1:]
            print(dtype)
            print(var_secondary.dimensions, type(var_secondary.dimensions))
    print(f'[INFO] Valid extra (secondary) extracts: {nvalid} of {ntotal}')


def main():
    path_main = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/Gustav_Dalen_Tower/IDEPIX/extracts'
    path_secondary = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/Gustav_Dalen_Tower/C2RCC/extracts'
    check(path_main, path_secondary)
    

if __name__ == '__main__':
    main()
