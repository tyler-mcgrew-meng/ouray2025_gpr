% directories in the “topo_corrections” project:

data: contains the GPR data axported from EkkoProject. The ‘csv’ subdirectory contains the data files exported from Ekko Project, which have been converted to the .mat files contained in the folder. These .mat files are composed of three matrix variables for each radargram: line#_t (M x 1), line#_x (1 X N), and line# (M x N)

figs: Folder to store output figures from matlab scripts

ouray_interps: Interpretation files exported from EkkoProject. The important variables contained here are the x and t coordinates for the interpretations in each radargram. 

output: Folder to write interpreted values with their respective coordinates

topo: Folder to store the topographic profiles and the corresponding E/N coordinates extracted with QGIS profile tool

trace_header: Potentially useful metadata for each trace collected in the survey
