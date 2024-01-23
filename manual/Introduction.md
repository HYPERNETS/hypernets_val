# **1. Introduction**

The Match-up Data Base (MDB) approach developed in the framework of the HYPERNETS project aimed at implementing validation analysis of satellite water products using HYPSTAR® data as reference.
 
The approach is based on the concept of MDB file that was first introduced by EUMETSAT, i.e., a NetCDF file including all the potential match-ups between satellite and in situ data within a time window. 

Moreover, it uses a set of open-source tools developed in Python to implement the validation analysis working with the MDB files. These MDB tools are divided into three modules: SAT_EXTRACT, MDB_builder and MDB_reader. All the modules are included in the hypernets_val repository available in GitHUB (https://github.com/HYPERNETS/hypernets_val).

This remainder of this document is organized as follows: 1) description of the MDB file structure (section 2); and 2) description of the MDB tools, with usage and examples (sections 3-5).

***

| [Table of contents](Index.md) | [2. MDB File Structure >](MDB_file_structure.md) |
| :------:| -----------:|
