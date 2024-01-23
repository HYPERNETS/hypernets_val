## ***2.3. MDB results (MDBr) files***

The MDB result (MDBr) files are generated using the GENERATEMU module of the MDB_reader module (section 5.1). In addition to the satellite and in situ data, it includes the match-ups for each specific wavelength after implementing the quality control. 

Dimensions and variables are summarized in Table 6 and Table 7. Global attributes are inherited from the MDB file (Table 3), with an update of the *creation_time* and *description* attributes.


Table 6: Dimensions of the MDBr files
|**Dimension**|**Description**|**Length**
|---|---|---
|*satellite_id*|Satellite measurements|Unlimited
|*satellite_bands*|Satellite bands|Depends on sensor/processor
|*rows*|y spatial coordinate|Default: 25
|*columns*|x spatial coordinate|Default: 25
|*insitu_id*|In situ measurements|Default: 40
|*insitu_original_bands*|In situ bands|1600 (for HYPSTAR®)
|*mu_id*|Match-up at wavelength|Unlimited

Match-up data require a new dimension (*mu_id*) defined as unlimited, with a current length equal to the number of satellite measurements (current length of *satellite_id*) by the number of satellite bands (maximum equal to the length of *satellite_bands*) included in the analysis. 


Table 7: Variables included in the MDBr files
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
|*time_difference*|Time difference|*satellite_id*





***

|[< 2.2 MDB files](MDB_files.md)| [Table of contents](Index.md) | [2.4 MDB concatenated results (MDBrc) files >](MDBrc_files.md) |
|:-----------| :------:| -----------:|