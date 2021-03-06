# Match-up DataBase for Validation of Satellite Data using HYPSTAR

The Match-up Data Base (MDB) includes match-ups between different satellite sensors and in situ data. It was developed as a way to improve the exchange, and processing of match-up data among different entities based on the concept by EUMETSAT (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html). There is one MDB file for each site and satellite sensor.

The workflow of the MDB creation:
- Create satellite extract
- Add _in situ_ data to the satellite extract
- Concatenation of all extract files to one single MDB file 

## Creation of Satellite Extract
The `hypernets_val/SAT_EXTRACT/sat_extract_OLCI.py` script creates the single satellite extracts in a NetCDF format. It can be run as:
```sh
cd SAT_EXTRACT
./sat_extract_OLCI.py -c CONFIG_FILE
```
where the config file has the parameters needed to select and extract the data from the satellite granules. An example of config file (hypernets_val/SAT_EXTRACT/config_sat_extract_OLCI.ini) used for the extract of EUMETSAT Sentinel-3/OLCI WFR data:
```
## config file for creating satellite extracts
[file_path]
sat_source_dir: /path/to/sat/data
output_dir: /path/where/to/save/output

[Time_and_sites_selection]
#first and last dates for analysis. Format: YYYY-MM-DD
time_start: 2021-01-01
time_stop: 2021-07-07
sites: VEIT

[satellite_options]
sat_prefix: satellite_
satellite: S3
platform: A,B
sensor: olci
# WRR o WFR for olci
resolution: WFR
# processor version. Note: without zeros to the left
proc_version: 7.00
extract_size: 25
window_size: 3
n_bands: 16
#apply BRDF to OLCI reflectance
BRDF: F
```

## Creation of the MDB file
The `hypernets_val//MDB_builder/MDB_builder.py` adds the _in situ_ data to the satellite extracts and concatenates them to a single MDB file in NetCDF format. It can be run as:
```sh
cd MDB_builder
MDB_builder.py -c CONFIG_FILE
```
where the config file has the parameters needed to select the satellite extracts and the in situ data. An example of config file (hypernets_val//MDB_builder/cconfig_MDB_builder_OLCI_HYPSTAR.ini.template) used for creating a MDB for EUMETSAT Sentinel-3/OLCI WFR data and HYPSTAR data:
```
## config file for MDBs builder
[file_path]
sat_extract_dir: /path/to/sat/extracts
ins_source_dir: /path/to/in/situ/data
output_dir: /path/where/to/save/output

[Time_and_sites_selection]
#first and last dates for analysis. Format: YYYY-MM-DD
time_start: 2021-01-01
time_stop: 2021-07-07
#insitu_type: PANTHYR, AERONET
insitu_type: HYPERNETS
#AERONET-OC  or PANTHYR sites:  ALL or any from AERONET-OC list. Default is ALL
sites: BEFR

[satellite_options]
sat_prefix: satellite_
satellite: S3
platform: A,B
sensor: olci
# WRR o WFR for olci
resolution: WFR
# processor version. Note: without zeros to the left
proc_version: 7.00
extract_size: 25
window_size: 3
n_bands: 16
#apply BRDF to OLCI reflectance
BRDF: F

```