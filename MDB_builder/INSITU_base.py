from netCDF4 import Dataset


class INSITUBASE:

    def __init__(self, mdb_options):
        self.mdb_options = mdb_options
        self.new_MDB = None

    def add_insitu(self, extract_path, ofile):
        self.start_add_insitu(extract_path, ofile)

    def start_add_insitu(self, extract_path, ofile):
        self.new_MDB = self.copy_nc(extract_path, ofile)

        # time_window = 2  # del mdb_options
        # self.new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds
        self.new_MDB.time_diff = self.mdb_options.insitu_options['time_max']

        n_insitu_id = self.mdb_options.insitu_options['n_insitu_id']
        n_insitu_bands = self.mdb_options.insitu_options['n_insitu_bands']
        self.new_MDB.createDimension('insitu_id', n_insitu_id)
        self.new_MDB.createDimension('insitu_original_bands',  n_insitu_bands)

        ##TIME VARIABLE
        insitu_time = self.new_MDB.createVariable('insitu_time', 'f4', ('satellite_id', 'insitu_id',), zlib=True,
                                                  complevel=6)
        insitu_time.units = "Seconds since 1970-1-1"
        insitu_time.description = 'In situ time in ISO 8601 format (UTC).'

        ##FILENAME VARIABLE
        insitu_filename = self.new_MDB.createVariable('insitu_filename', 'S2', ('satellite_id', 'insitu_id'), zlib=True,
                                                      complevel=6)
        insitu_filename.description = 'In situ filename.'

        #ORIGINAL BANDS VARIABLE
        insitu_original_bands = self.new_MDB.createVariable('insitu_original_bands', 'f4', ('insitu_original_bands'),
                                                            fill_value=-999, zlib=True, complevel=6)
        insitu_original_bands.description = 'In situ bands in nm.'

        #RRS VARIABLE
        insitu_Rrs = self.new_MDB.createVariable('insitu_Rrs', 'f4', ('satellite_id', 'insitu_original_bands', 'insitu_id'),
                                            fill_value=-999, zlib=True, complevel=6)
        insitu_Rrs.description = 'In situ Rrs'

        #TIME DIFFERENCE VARIABLE
        time_difference = self.new_MDB.createVariable('time_difference', 'f4', ('satellite_id', 'insitu_id'),
                                                 fill_value=-999,
                                                 zlib=True, complevel=6)
        time_difference.long_name = "Absolute time difference between satellite acquisition and in situ acquisition"
        time_difference.units = "seconds"

        #self.new_MDB.close()

    def copy_nc(self, ifile, ofile):
        with Dataset(ifile) as src:
            dst = Dataset(ofile, 'w', format='NETCDF4')
            # copy global attributes all at once via dictionary
            dst.setncatts(src.__dict__)
            # copy dimensions
            for name, dimension in src.dimensions.items():
                dst.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else None))
            # copy all file data except for the excluded
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                # copy variable attributes all at once via dictionary
                dst[name].setncatts(src[name].__dict__)

                dst[name][:] = src[name][:]
        return dst
