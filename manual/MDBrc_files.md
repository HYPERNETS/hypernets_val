## ***2.4. MDB concatenated results (MDBrc) files***

The MDB concatenated results (MDBrc) are generated using the mode CONCATENATE of the MDB_reader module (section 5.2.). They combine results from different MDBr files, adding a set of *flag_* variables to identify each match-up.  It uses the same dimensions as MDBr files (Table 6). Variables are summarized in Table 8. 

Table 8: Variables included in the MDBrc files
|**Variable**|**Description**|**Dimensions**
|---|---|---
|*satellite_bands*|Band wavelengths (nm)|*satellite_bands*
|*satellite_time*|Overpass time|*satellite_id*
|*satellite_Rrs*|Satellite-derived *Rrs*|*satellite_id, satellite_bands, rows, columns*
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
|*insitu_Rrs*|In situ *Rrs*|*satellite_id, insitu_original_bands, insitu_id*
|*insitu_Rrs_nosc*|In situ *Rrs* without correction for the NIR similarity spectrum|*satellite_id, insitu_original_bands, insitu_id*
|*insitu_quality_flag*|Quality Flag Dataset|*satellite_id, insitu_id*
|*insitu_site_flag*|Site Flag Dataset|*satellite_id, insitu_id*
|*insitu_OAA*|Observation Azimuth Angle|*satellite_id, insitu_id*
|*insitu_OZA*|Observation Zenith Angle|*satellite_id, insitu_id*
|*insitu_SAA*|Sun Azimuth Angle|*satellite_id, insitu_id*
|*insitu_SZA*|Sun Zenith Angle|*satellite_id, insitu_id*
|*mu_ins_rrs*|Match-up in situ *Rrs*|*mu_id*
|*mu_sat_rrs*|Match-up satellite *Rrs*|*mu_id*
|*mu_wavelength*|Match-up wavelength|*mu_id*
|*mu_satellite_id*|Match-up *satellite_id*|*mu_id*
|*mu_valid*|Match-up validity|*satellite_id*
|*mu_insitu_id*|Match-up *insitu_id*|*satellite_id*
|*mu_ins_time*|Match-up in situ time|*satellite_id*
|*mu_sat_time*|Match-up satellite time|*satellite_id*
|*mu_time_diff*|Match-up time difference |*satellite_id*
|*flag_ac*|Atmospheric correction|*satellite_id*
|*flag_site*|Site|*satellite_id*
|*flag_satellite*|Satellite mission|*satellite_id*
|*flag_sensor*|Satellite sensor|*satellite_id*
|*time_difference*|Time difference|*satellite_id*

The new flag variables use the following flags included in the global attributes of the MDBr files:
*flag_ac: satellite_aco_processor*
*flag_site: insitu_site_name*
*flag_satellite: satellite + platform*
*flag_sensor: sensor*

A value is assigned to each flag using 2<sup>n</sup> with *n* being consecutive numbers staring from 0 (i.e. flag values would be 1, 2, 4, 8, etc). 

Global attributes are inherited from the MDBr files (Table 3). In case of attributes used for flagging (*satellite_aco_processor*, *insitu_site_name*, *satellite*, *platform*, *sensor*), as well as *insitu_lat* and *insitu_lon*, attributed values are updated using a list (comma-separated values) with all the values inherited of the MRBr files included in the concatenation. *creation_time* and *description* attributes are also updated. 




***

|[< 2.3 MDB results (MDBr) files](MDBr_files.md)| [Table of contents](Index.md) | [3. SAT_extract >](SAT_extract.md) |
|:-----------| :------:| -----------:|