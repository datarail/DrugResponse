function Design = TextDesignFile_importer(filename)
% Design = TextDesignFile_importer(filename)
%   Read a tsv file and convert the information in an array of Design structure
%   with standard fields for treatment. Assume only one plate therefore a
%   single entry per well.
%
%   filename :  path and file name to a treatment saved as tsv file with
%               the following columns:
%                   - DrugName (and HMSLid)
%                   - Well
%                   - Conc (for concentration in uM)
%               and annotation columns:
%                   - DMSO (added to the treated_wells)
%               and optional columns (stored as 'Perturbations')
%                   - pert_type (e.g. trt_cp, vehicle_ctl, ...)
%                   - SeedingNumber
%                   - other perturbations (e.g. EGF/...)
%                   - DrugName2/3/... and Conc2/3/... for multiple drugs per well
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO)
%                   - well_volume (in uL)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)
%


assert(exist(filename, 'file')>0, 'Design file %s missing', filename)

t_design = tsv2table(filename);

Design = TreatmentTableToDrugDesign(t_design);
