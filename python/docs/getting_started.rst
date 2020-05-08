Getting started
===============

.. |dissapointed| replace:: 😞

.. |eyes| replace::


Working with data generated from Columbus
-----------------------------------------

 - Object level data for a given well is saved as a ``.txt`` file. The data for all wells in a given plate is saved in a folder. The name of the folder typically takes the form ``abcdef[123]``. The ``abcdef`` prefix in the folder name should correspond to the barcode assigned to the plate.
 - ``cd`` into the directory that contains the object level data folders.
 -  Start a Jupyter notebook or Ipython session and execute the lines of code below:

::

   from cell_cycle_gating import run_cell_cycle_gating as rccg
	
   obj_folder = 'abcdef[123]' # Name of object level folder
	
   # Map user defined channel names to standarized names required by the script
   ndict = {'Nuclei Selected - EdUINT': 'edu',
            'Nuclei Selected - DNAcontent': 'dna',
	    'Nuclei Selected - LDRTXT SER Spot 8 px' : 'ldr',
	    'Nuclei Selected - pH3INT': 'ph3'}

   # Run gating script	
   dfs = rccg.run(obj_folder, ndict)

- The dataframe ``dfs`` returns well-level summary of number of live/dead cells and fraction of cells in each phase of the cell cycle.
- The script saves a pdf showing the gating on each DNA v EDU scatter plot for review. By default the pdf file uses the name of the folder as the file name `i.e` ``summary_abcdef[123].pdf``
- The dataframe df is also saved as a .csv file with the same name as the object level folder `i.e` ``summary_abdef[123].csv``
  

Working with data generated from IXM
------------------------------------

- The standard column names generated in  Metaexpress differ from the Operetta. Please update ``ndict`` as below and add an additional argument to ``rccg.run()``. Further, the entire dataset for a given plate is saved as a single .txt file instead of seperate files per well in a folder. See below:

  
::

   from cell_cycle_gating import run_cell_cycle_gating as rccg
	
   obj_file = 'filename.txt' # Name of object level file
   

   # Map user defined channel names to standarized names required by the script
   ndict = {'Well Name' : 'well',
            'Cell: EdUrawINT (DDD-bckgrnd)' : 'edu',
            'Cell: LDRrawINT (DDD-bckgrnd)' : 'ldr',
            'Cell: DNAcontent (DDD-bckgrnd)' : 'dna'}

   # Run gating script	
   dfs = rccg.run(obj_file, ndict, system='ixm', header=7)


Merging metadata information to output
--------------------------------------

If you have well level metadata that maps each well to sample conditions, then the above code block can be modified as follows:

::

   # Load metadata file
   import pandas as pd
   dfm = pd.read_csv('metadata.csv')

   # Run gating script, this time passing dfm as an additional argument.
   dfs = rccg.run(obj, ndict, dfm)

Note that the metadata file should contain the following header columns:
   - ``barcode``, ``well``, ``cell_line``, ``agent``, ``concentration``.
   - Fields in the ``barcode`` column should match the prefix in the folder name. `i.e` ``abdcdef``


Not really **Deep** Dye Drop
----------------------------
By default, the gating code expects that you have all 4 channels `i.e` DNA, EdU, LDR, pH3. However, if you do not have LDR and/or pH3 channels, modify the main line of the code as shown below:

::
     
   dfs = rccg.run(obj, ndict, dfm,
                  ph3_channel=False, # If no pH3 channel
		  ldr_channel=False # If no LDR channel
		 )
		

Gating corrections using control wells
--------------------------------------

Automated gating does not always work well, in which case you can apply the automated gating from control wells for a given cell line across corresponding treatment wells. In the first line of code below, by setting ``control_based_gating=True``, automated gating is run only on the control plates. The results are saved in .csv and .pdf files with the prefix ``control_summary_``. The function also returns a second dataframe ``dfg`` that contains information on the gates in the control wells. In the second line, the script is run a second time with the arguement ``control_gates=dfg`` so that errors in automated gating are corrected based on control gating.

::

   dfs, dfg = rccg.run(obj, ndict, dfm,
                       control_based_gating=True)
   dfs2 = rccg.run(obj, ndict, dfm, control_gates=dfg)

If you want to manually adjust gates across all wells, you can provide a list of fudge factors `i.e.` by what magnitude and in which direction you want to change the DNA gates. There are 4 gates you can adjust; G1-left, G1-right, G2-left, and G2-right. For instance, if you want to move G2-left (3rd gate) furrther left by a magnitude of 0.05, set ``fudge_gates=[0, 0, -0.05, 0]``. If you want to move G2-right (4th gate) by a magnitue of 0.2 to the right, set ``fudge_gates=[0, 0, 0, 0.2]``. In the code example below, we have applied gates from the control but also decided to move the 3rd and 4th gates to the left by 0.05 and 0.1 units in log(DNA) scale.

::

   dfs2 = rccg.run(obj, ndict, dfm, control_gates=dfg,
                   fudge_gates=[0, 0, -0.05, -0.1])
		  

		  
Additional plotting functions
-----------------------------
.. - To plot DNA content distributions:
..
..
..  import pandas as pd
..   from cell_cycle_gating import plot_dna_distributions

..   dfm = pd.read_csv('metadata.csv')
..   obj = 'abcdef[123]'

..   plot_dna_distributions(obj, dfm)

- To plot cell cycle fractions:
::

   import pandas as pd
   from cell_cycle_gating import plot_fractions

   dfs = pd.read_csv('summary_abcdef[123].csv')
   plot_fractions(dfs)
   
   
   
		       
