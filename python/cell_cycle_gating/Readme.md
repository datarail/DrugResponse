# Deep Dye Drop based cell cycle gating

cd into the directory that contains the object level data folders. Start a Jupyter notebook or Ipython session.
For each plate (object level folder), use the template below to compute live/dead status and cell cycle phases of individual cells.
The dataframe `df` returns well-level summary of number of live/dead cells and fraction of cells in each phase of the cell cycle.
``` 
from cell_cycle_gating import run_cell_cycle_gating as rccg
obj = 'path_to_object_level_data_folder'
df = rccg.run(obj)
```    

![Alt text](example_plots/example_plot.jpg?raw=true "Title")