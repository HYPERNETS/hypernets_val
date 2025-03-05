import os.path
import sys
main_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(main_path)

path_aceasy = os.path.join(main_path,'aceasy')
if os.path.exists(path_aceasy):
    print(f'[INFO] __init__: {path_aceasy} added to the path')
    sys.path.append(path_aceasy)