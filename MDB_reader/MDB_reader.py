from MDBFile import MDBFile
from MDBFileList import MDBFileList
import os


class MDB_READER():
    def __init__(self, path_mdb):
        if path_mdb is not None:
            self.mfile = MDBFile(path_mdb)


def main():
    path_base = '/mnt/c/DATA_LUIS/OCTAC_WORK/BAL_EVOLUTION/EXAMPLES/TRIMMED/MDBs'
    #name_mdb = 'MDB_S3A_OLCI_WFR_STANDARD_L2_AERONET_Gustav_Dalen_Tower.nc'
    name_mdb = 'MDB_S3A_OLCI_EFR_POLYMER_L2_AERONET_Gustav_Dalen_Tower.nc'
    path_mdb = os.path.join(path_base, name_mdb)
    reader = MDB_READER(path_mdb)
    reader.mfile.qc_sat.compute_invalid_masks(0)
    #reader.mfile.qc_sat.compute_flag_stats(5)
    #reader.mfile.qc_sat.compute_flag_masks(5)
    #reader.mfile.qc_sat.add_theshold_mask(0, -1, 0.001, 'greater')
    # reader.mfile.qc_sat.add_band_statistics(-1, 560,'CV', 20, 'greater')
    # b = reader.mfile.qc_sat.compute_statistics(0)
    # print(b)
    # b = reader.mfile.qc_sat.do_check_statistics()
    # print(b)




if __name__ == '__main__':
    main()
