## ***3. SAT_extract***

The satellite data are provided as satellite extracts created using the tools available in the SAT_EXTRACT module. Different Python extracts tools were created for working with different satellite sensors and/or processors (Table 9).

Table 9: Satellite extract tools implemented in the *hypernets_val* repository
|**SAT_EXTRACT tool**|**Satellite**|**Processor**
|---|---|---
|sat_extract_OLCI|S3A, S3B|STANDARD (WFR)
|sat_extract_CMEMS|CMEMS L3 and L4 products|CMEMS
|sat_extract_POLYMER|S3A, S3B, S2A, S2B|POLYMER
|sat_extract_ACOLITE|S3A, S3B, S2A, S2B|ACOLITE
|sat_extract_C2RCC|S3A, S3B, S2A, S2B|C2RCC
|sat_extract_NASA|MODIS, VIIRS, VIIRSJ|STANDARD (L2)

The format of the required source files is NetCDF for all the extracts tools except for sat_extract_OLCI, which is based on the Sentinel-3 SAFE format.

All the extract tools are run using a script passing as argument a configuration file with all the parameters and options:

```
python sat_extract_OLCI.py -c extract_config.ini -v
```
The tool creates an extract file for each available source image covering the specified in situ location within a temporal range, which can be defined using the start and stop dates or with a date list. Source files can be optionally filtered using a wild card expression. 

Output extract files are NetCDF files containing satellite data from a specific product and for an extraction window (by default 25 x 25 pixels) centred on the specified in situ location. It includes a set of attributes defining the site (*insitu_site_name, insitu_lat, insitu_lon*) and satellite product (*satellite, platform, sensor, resolution, satellite_aco_processor, satellite_proc_version*). 

Although satellite attributes are retrieved from the satellite source files, values (except for *satellite_aco_processor*) can optionally be established using the configuration file. 

Configurations files are organized in three sections: **file_path**, **Time_and_sites_selection** and **satellite_options**. 

**[file_path]**

***sat_source_dir:*** Path to the directory including the satellite source files with a specific format depending on the extract tool.  Required. 

***sat_source_dir_organization:*** Structure of the source directory in case of files organized in sub-folders indicating the date. It uses YYYY for year, mm for month, dd for day of the month and jjj for the Julian date. For instance, YYYY/jjj or YYYY/mm/dd. Optional. 

***output_dir:*** Output folder for the satellite extract files. Required.

***tmp_dir:*** Temporary folder to decompress source files in compressed formats (i.e., zip or tar).  Decompressed files are deleted after creating the extract. Optional (*sat_source_dir* is used to decompress source files if this option is not available). 

**[Time_and_sites_selection]**

***time_start:*** First date for analysis. Format: YYYY-mm-dd. It used in combination with *time_stop*.

***time_stop:*** Last date for analysis. Format: YYYY-mm-dd. It used in combination with *time_start*. 

***time_list_file:*** Path to a text file including a date list in format YYYY-mm-dd. If this option is given, *time_start* and *time_stop* parameters are not used. 

***site:*** Site name. The coordinates of the following WATERHYPERNET sites are already included in the tool: VEIT, GAIT, BEFR, MAFR, M1BE, LPAR. For other sites, please provide latitude and longitude using *insitu_lat* and *insitu_lon*. 

***insitu_lat:*** Latitude of the site (required if site is not a WATERHYPERNET site)

***insitu_lon:*** Longitude of the site (required if site is not a WATERHYPERNET site)

**[satellite_options]**

***extract_size:*** Size of the extraction window. Optional. Default: 25 pixels. 

***wce:*** Wild card expression to filter the source files based on their file name. For instance, S3A* limits the extracts to Sentinel-3A. Optional. Default: None.

***BRDF:*** If True, it applies BRDF correction (only for sat_extract_OLCI). Boolean. Default: False. Optional. 

***satellite:*** Satellite attribute. Optional. 

***platform:*** Platform attribute. Optional. 

***sensor:*** Sensor attribute. Optional.

***resolution:*** Resolution attribute. Optional. 

***satellite_proc_version:*** Atmospheric correction version. Optional. 

Example of satellite extract configuration file: 

~~~
## config file for creating satellite extracts
[file_path]
sat_source_dir: /store3/SAT_EXTRACTS/OLCI/source
output_dir: /store3/SAT_EXTRACTS/OLCI/extracts
tmp_dir: /store3/SAT_EXTRACTS/tmp
sat_source_dir_organization: YYYY/jjj
 
[Time_and_sites_selection]
time_start: 2019-01-01
time_stop: 2023-03-31
site: BEFR
time_list_file: /store3/SAT_EXTRACTS/configFiles/date_lists/list_BEFR.txt

[satellite_options]
extract_size: 25
BRDF: F
wce:  S3A*
~~~

***

|[< 2.4 MDB concatenated results (MDBrc) files](MDBrc_files.md)| [Table of contents](Index.md) | [4. MDB_builder >](MDB_builder.md) |
|:-----------| :------:| -----------:|