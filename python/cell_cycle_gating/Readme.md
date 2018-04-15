# Deep Dye Drop based cell cycle gating

## Requirements
* peakutils (available on [PyPI](https://pypi.python.org/pypi/PeakUtils)). It can be installed from command line as shown below
```
$ pip install peakutils
```
## Installation
* The cell cycle gating repository can be installed from command line as shown belo
```
$ git clone https://github.com/datarail/DrugResponse.git
```
* add `datarail/DrugResponse/python/cell_cycle_gating` to your `PYTHONPATH` to enable importing the modules under `cell_cycle_gating` from any location on your machine

## Getting started
* cd into the directory that contains the object level data folders. 
* Start a Jupyter notebook or Ipython session.
* For each plate (object level folder), use the template below to compute live/dead status and cell cycle phases of individual cells.
* The dataframe `df` returns well-level summary of number of live/dead cells and fraction of cells in each phase of the cell cycle.
* The script saves a pdf showing the gating on each DNA v EDU scatter plot for review. 
* The dataframe `df` is also saved as a .csv file with the same name as the object level folder.
``` 
from cell_cycle_gating import run_cell_cycle_gating as rccg
obj = 'path_to_object_level_data_folder'
df = rccg.run(obj)
```    

![Alt text](example_plots/example_plot.jpg?raw=true "Title")
