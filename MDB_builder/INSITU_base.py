from netCDF4 import Dataset


class INSITUBASE:

    def __init__(self,mdb_options):
        self.mdb_options = mdb_options
        self.new_mdb = None

    def add_insitu(self,extract_path,ofile):
        self.start_add_insitu(extract_path,ofile)

    def start_add_insitu(self,extract_path,ofile):
        self.new_MDB = self.copy_nc(extract_path, ofile)

        time_window = 2  #del mdb_options
        self.new_MDB.time_diff = f'{time_window * 60 * 60}'  # in seconds

        self.new_MDB.createDimension('insitu_id', 30)
        self.new_MDB.createDimension('insitu_original_bands', 1605)

        self.new_MDB.close()

    def copy_nc(ifile, ofile):
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