# **2. MDB file structure**

MDB files are based on the open-source NetCDF-4 (Network Common Data Form) file format (https://www.unidata.ucar.edu/software/netcdf), which is built on top of the Hierarchical Data Format version 5 (HDF5) (https://hdfgroup.github.io/hdf5). It supports large, complex and heterogenous data by using a directory-like structure to organize the data within the file. Data in MDB files are stored in variables, i.e., multi-dimensional arrays of values of the same type.

MDB files are built in different stages, adding new data:

a. Satellite extract files. They are generated using the SAT_extract module. They only include satellite data (section 3).

b. MDB files: They are produced using the MDB_builder module. The include satellite and in situ data (section 4).

c. MDB result files (MDBr): They are produced using the MBD_reader module (mode GENERATEMU) and include the match-ups between the satellite and situ data after implementing the quality control (section 5.1). 

d. MDB concatenated files (MDBrc): They are generated using the MDB_reader module (mode CONCATENATE). The put together MDB results files adding a set of flag variables. 

***

|[< 1. Introduction](Introduction.md)| [Table of contents](Index.md) | [2.1 Satellite extract files >](sat_extract_structure.md) |
|:-----------| :------:| -----------:|