import os.path
import sys
main_path = os.path.dirname(os.path.dirname(__file__))
sys.path.append(main_path)
path_mdb_reader = os.path.join(main_path,'MDB_reader')
sys.path.append(path_mdb_reader)
