## ***2.2. MDB files***

MDB files are generated using the MDB_builder module (section 4) combining satellite data from the extract files (section 3) and in situ data from HYPSTAR® L2 files for the same site and within a time window. Satellite data are directly inherited from the satellite extract files (section 2.1), whereas two new dimensions are required to describe the in situ data (*insitu_id* and *insitu_original_bands*). Note that the current length of *satellite_id* dimension (defined as unlimited) is the number of satellite measurements (i.e. satellite extracts) included in the MDB file, whereas the length of *insitu_id* dimension would be equal to the maximum number of in situ measurements (*Rrs* spectra) that could be associated with a specific satellite measurement (by default is set to 40).

Dimensions and variables are summarized in Table 4 and Table 5. Global attributes are inherited from the satellite extract files (Table 3), with an update of the *creation_time* and *description* attributes. 

Table 4: Dimensions of the MDB files
|**Dimension**|**Description**|**Length**
|---|---|---
|*satellite_id*|Satellite measurements|Unlimited
|*satellite_bands*|Satellite bands|Depends on sensor/processor
|*rows*|y spatial coordinate|Default: 25
|*columns*|x spatial coordinate|Default: 25
|*insitu_id*|In situ measurements|Default: 40
|*insitu_original_bands*|In situ bands|1600 (for HYPSTAR®)


Table 5: Variables included in the MDB files
|**Variable**|**Description**|**Dimensions**
|---|---|---
|*satellite_bands*|Band wavelengths (nm)|*satellite_bands*
|*satellite_time*|Overpass time|*satellite_id*
|*satellite_Rrs*|Satellite-derived Rrs|*satellite_id, satellite_bands, rows, columns*
|*satellite_latitude*|Latitude|*satellite_id, rows, columns*
|*satellite_longitude*|Longitude|*satellite_id, rows, columns*
|*satellite_AOT_0865p50*|Aerosol Optical Thickness|*satellite_id, rows, columns*
|*satellite_flag*|Flags Data Set|*satellite_id, rows, columns* 
|*satellite_OAA*|Observation Azimuth Angle|*satellite_id, rows, columns*
|*satellite_OZA*|Observation Zenith Angle|*satellite_id, rows, columns*
|*satellite_SAA*|Sun Azimuth Angle|*satellite_id, rows, columns*
|*satellite_SZA*|Sun Zenith Angle|*satellite_id, rows, columns*
|*insitu_original_bands*|Instrument wavelenghts (nm)|*insitu_original_bands*
|*insitu_time*|Measurement time|*satellite_id, insitu_id*
|*insitu_Rrs*|In situ Rrs|*satellite_id, insitu_original_bands, insitu_id*
|*insitu_Rrs_nosc*|In situ Rrs without correction for the NIR similarity spectrum|*satellite_id, insitu_original_bands, insitu_id*
|*insitu_quality_flag*|Quality Flag Dataset|*satellite_id, insitu_id*
|*insitu_site_flag*|Site Flag Dataset|*satellite_id, insitu_id*
|*insitu_OAA*|Observation Azimuth Angle|*satellite_id, insitu_id*
|*insitu_OZA*|Observation Zenith Angle|*satellite_id, insitu_id*
|*insitu_SAA*|Sun Azimuth Angle|*satellite_id, insitu_id*
|*insitu_SZA*|Sun Zenith Angle|*satellite_id, insitu_id*
|*time_difference*|Time difference|*satellite_id*





***

|[< 2.1 Satellite extract files](sat_extract_structure.md)| [Table of contents](Index.md) | [2.3 MDB results (MDBr) files >](MDBr_files.md) |
|:-----------| :------:| -----------:|