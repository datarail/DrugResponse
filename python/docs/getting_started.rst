Getting started
===============

.. |dissapointed| replace:: 😞

.. |eyes| replace:: https://www.smileysapp.com/emojis/watching-you.png


Quickstart
----------

 - Object level data for a given well is saved as a ``.txt`` file. The data for all wells in a given plate is saved in a folder. The name of the folder typically takes the form ``abcdef[123]``. The ``abcdef`` prefix in the folder name should correspond to the barcode assigned to the plate.
 - ``cd`` into the directory that contains the object level data folders.
 -  Start a Jupyter notebook or Ipython session and execute the lines of code below:

::

   from cell_cycle_gating import run_cell_cycle_gating as rccg
	
   obj = 'abcdef[123]' # Name of object level folder
	
   # Map user defined channel names to standarized names required by the script
   ndict = {'Nuclei Selected - EdUINT': 'edu',
            'Nuclei Selected - DNAcontent': 'dna',
	    'Nuclei Selected - LDRTXT SER Spot 8 px' : 'ldr',
	    'Nuclei Selected - pH3INT': 'ph3'}

   # Run gating script	
   dfs = rccg.run(obj, ndict)

- The dataframe ``dfs`` returns well-level summary of number of live/dead cells and fraction of cells in each phase of the cell cycle.
- The script saves a pdf showing the gating on each DNA v EDU scatter plot for review. By default the pdf file uses the name of the folder as the file name `i.e` ``summary_abcdef[123].pdf``
- The dataframe df is also saved as a .csv file with the same name as the object level folder `i.e` ``summary_abdef[123].csv``


Merging metadata information to output
--------------------------------------

If you have (and you should have |eyes|) well level metadata that maps each well to sample conditions, then the above code block can be modified as follows:

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
		

Gating only on controls |dissapointed|
--------------------------------------

Automated gating does not always work well, in which case you can apply the automated gating from control wells for a given cell line across corresponding treatment wells.

::

   dfs = rccg.run(obj, ndict, dfm,
                  control_based_gating=True)

		  
Additional plotting functions
-----------------------------
By default, the script outputs individual scatter plots for each well in a plate.

 - To plot DNA content distributions:
::

   import pandas as pd
   from cell_cycle_gating import plot_dna_distributions

   dfm = pd.read_csv('metadata.csv')
   obj = 'abcdef[123]'

   plot_dna_distributions(obj, dfm)

- To plot cell cycle fractions:
::

   import pandas as pd
   from cell_cycle_gating import plot_fractions

   dfs = pd.read_csv('summary_abcdef[123].csv')
   plot_fractions(dfs)
   
   
   
		       
