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
python sat_extract_OLCI.py -c *extract_config.ini* -v
```

***

|[< 2.4 MDB concatenated results (MDBrc) files](MDBrc_files.md)| [Table of contents](Index.md) | [4. MDB_builder >](MDB_builder.md) |
|:-----------| :------:| -----------:|