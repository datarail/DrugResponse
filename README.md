# DrugResponse
Analysis of drug response based on cell staining

## Deep Dye Drop based cell cycle gating

### Installation via conda/mamba
If necessary, install mamba:
```
conda install mamba
```
Create a conda environment with the necessary dependencies:
```
mamba env create -n deep_dye_drop --file=env_conda_2025.yaml
```
Activate the new conda environment
```
conda activate deep_dye_drop
```

Open Jupyter Lab:
```
jupyter-lab
```

If Jupyter Lab doesn't automatically open in your browser, visit http://localhost:8888/

Select the "deep-dye-drop" conda environment for the Python kernel.

### Installation via docker

Pull the docker image
```
docker pull labsyspharm/deep-dye-drop
```

Change to the directory you would like to mount via docker
```
cd <your directory>
```

Create a docker container running on port 7777 (you may choose any free port)
```
docker run -d -p 7777:8888 -v "${PWD}":/home/jovyan/work --name ddd_notebook labsyspharm/deep-dye-drop start-notebook.sh --IdentityProvider.token=''
```
 This will mount the current working directory to the "work" folder within the container, giving you access to the files and jupyter notebooks therein.

 Note: "jovyan" is the default username used by the [Jupyter Docker Stacks](https://jupyter-docker-stacks.readthedocs.io/en/latest/) docker images upon which our docker image is built.

 Open Jupyter Lab in your browser at http://localhost:7777/ or on the port of your choosing.

 Select the "deep-dye-drop" conda environment for the Python kernel.

### Cell cycle gating example

See the "python/cell_cycle_gating/examples/DDR_example.ipynb" notebook.

### Getting started

* cd into the directory that contains the object level data folders. 
* Start a Jupyter notebook or Ipython session.
* For each plate (object level folder), use the template below to compute live/dead status and cell cycle phases of individual cells.
* The dataframe `df` returns well-level summary of number of live/dead cells and fraction of cells in each phase of the cell cycle.
* The script saves a pdf showing the gating on each DNA v EDU scatter plot for review. 
* The dataframe `df` is also saved as a .csv file with the same name as the object level folder.
  ```python 
  from cell_cycle_gating import run_cell_cycle_gating as rccg
  obj = 'path_to_object_level_data_folder'
  df = rccg.run(obj)
  ```    

![Alt text](python/cell_cycle_gating/example_plots/example_plot.jpg?raw=true "Title")

### LICENSE and FUNDING
Deep Dye Drop's automated gating package is currenlty available under the MIT license .The package was developed with funding from U54 grant HL127365, "The Library of Integrated Network-Based Cellular Signatures" under the NIH Common Fund program, and NCI U54 grant CA225088 for the Harvard Medical School (HMS) Center for Cancer Systems Pharmacology (CCSP).
