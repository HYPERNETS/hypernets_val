# Match-up DataBase for Validation of Satellite Data using HYPSTAR_

The Match-up Data Base (MDB) includes match-ups between different satellite sensors and in situ data. It was developed as a way to improve the exchange, and processing of match-up data among different entities based on the concept by EUMETSAT (https://ocdb.readthedocs.io/en/latest/ocdb-MDB-user-manual.html). There is one MDB file for each site and satellite sensor.

The workflow of the MDB creation:
- Create satellite extract
- Add _in situ_ data to the satellite extract
- Concatenation of all extract files to one single MDB file 

## Creation of Satellite Extract
The `sat_extract_OLCI.py` script creates the single satellite extracts in a NetCDF format. It can be run as:
``` 
cd SAT_EXTRACT
./sat_extract_OLCI.py -c CONFIG_FILE
```
where the config file has the parameters needed to select and extract the data from the satellite granules.

