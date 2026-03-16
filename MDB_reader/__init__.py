import os,sys
hypernets_path = os.path.dirname(os.path.dirname(__file__))
print(f'[INFO] __init__: Adding hypernets path: {hypernets_path}')
sys.path.append(hypernets_path)
main_path = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
print(f'[INFO] __init__: Adding main path: {main_path}')
sys.path.append(main_path)

path_aceasy = os.path.join(main_path,'aceasy')
if os.path.exists(path_aceasy):
    print(f'[INFO] __init__: {path_aceasy} added to the path')
    sys.path.append(path_aceasy)