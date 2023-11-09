# d3v-sgd
Design visualizer for  Ship Grillage Design

## Instalation of environment 
Run miniconda or anaconda command prompt/power shell  
conda create -n d3vsgd_env python=3.9 -c conda-forge  
conda activate d3vsgd_env  
pip install pyside6  
pip install openmesh  
pip install pyopengl  
pip install scipy  

clone d3v repository: (e.g. https://github.com/linaetal-fsb/d3v.git)  
clone d3v-sgd repository: (e.g. https://github.com/pprebeg/d3v-sgd.git)  

## Usage
Run as: Path_to_d3v\src\main.py -a Path_to_d3v-sgd<br>
e.g. if d3v is located in D:/dev/d3v, and d3v-sgd in D:/dev/d3v-sgd: <br>
python D:/dev/d3v/main.py -a D:/dev/d3v-sgd

Note: 1 "-a Path_to_d3v-sgd" is the argument
2 the "commands" folder needs to be the subfolder in d3v-sgd folder

## List of contributors:
Gordan Kos
Tomislav Pavlović
Pero Prebeg
