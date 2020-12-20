import pandas as pd
import sys
from cell_cycle_gating import run_cell_cycle_gating as rccg

metadata = sys.argv[1]
plate = sys.argv[2]

ndict = {'Nuclei Selected - DNAcontent': 'dna',
         'Nuclei Selected - LDRINT' : 'ldrint',
         'WellName': 'well',
        'Nuclei Selected - EdUINT': 'edu',
         'Nuclei Selected - pH3INT': 'ph3'}

#ndict = {'Cell: DNAcontent (DDD-bckgrnd)' : 'dna',
#         'Cell: EdUrawINT (DDD-bckgrnd)': 'edu',
#         'Well Name': 'well'}



dfm = pd.read_csv(metadata)
dfm.loc[dfm.role == 'negative_control', 'agent'] = 'DMSO'
dfm['timepoint'] = dfm['timepoint'].replace([0], 'time0_ctrl')
time0_plates = dfm[dfm.timepoint == 'time0_ctrl'].barcode.unique()

barcode = plate.split('.txt')[0].split('[')[0]
#barcode = plate.split('.txt')[0]

# Run control gating for treatment plates only
if barcode in time0_plates:
    dfs = rccg.run(plate, ndict, dfm=dfm, fudge_gates=[0, 0.01, -0.01, 0])
else:
    dfs, dfg = rccg.run(plate, ndict, dfm=dfm, control_based_gating=True)
    dfs = rccg.run(plate, ndict, dfm=dfm, control_gates=dfg, fudge_gates=[0, 0.01, -0.01, 0])
