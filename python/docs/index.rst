.. Deep Dye Drop gating documentation master file, created by
   sphinx-quickstart on Thu Dec  6 10:17:51 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for automated cell cycle gating
=============================================

The Deep Dye Drop assay provides detailed readouts on the cell cycle in addition to the live/dead counts.
Cells are treated with Hoechst (to stain nuclei) and LIVE/DEAD Red (LDR, to stain dead cells) dyes in multi-well plates.
In addition, cells can be treated with EdU (to stain S-phase cells) and an antibody for phospho-histone H3 (to stain M-phase cells).
Live cells are assigned to a phase of the cell cycle based on their DNA content, and EdU and pH3 intensities. Subsequently, cells are
imaged at the single-cell level. The `cell_cycle_gating` script described in this documentation takes the resultant single cell imaging
data as input and automatically classifies cells into live or dead and further classifies live cells into G1, G2, S, or M phase of the cell cycle.


.. toctree::
   :maxdepth: 3

   license.rst
   installation.rst
   getting_started.rst
   modules/index  


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
