function hpdd_exporter(hpdd_filename, Designs, t_plateinfo)
%hpdd_exporter(hpdd_filename, Designs, t_plateinfo)
%   Write an hpdd file based on array of Designs and a plate barcode table.
%
%   base_pathname : path and file name to a D300 protocol file ('.hpdd'
%                       automatically appended -- do not include it here)
%   Designs :       array of design structures with the following fields:
%                       - plate_dims (plate dimension)
%                       - treated_wells (wells treated to be backfilled)
%                           by default the backfill is the same as Drugs
%                           (or DMSO), but it can specified as:
%                               0: DMSO
%                               1: Aqueous + Brij35
%                               2: Aqueous + Brij35 + glycerol
%                               3: Aqueous + Tryton X100
%                               4: Aqueous + Tryton X100 + glycerol
%                               5: Aqueous + Tween 20
%                               6: Aqueous + Tween 20 + glycerol
%                       - well_volume (in uL)
%                       - Drugs (structure with DrugName, HMSLid, stock_conc
%                           and layout - concentration given in uM)
%                           Vehicle of the drug can be specified as an
%                           additional field 'Vehicle' with values:
%                               0: DMSO (default)
%                               1: Aqueous + Brij35
%                               2: Aqueous + Brij35 + glycerol
%                               3: Aqueous + Tryton X100
%                               4: Aqueous + Tryton X100 + glycerol
%                               5: Aqueous + Tween 20
%                               6: Aqueous + Tween 20 + glycerol
%
%   t_plateinfo :   table (or file name to a tsv table) with columns:
%                       - Barcode
%                       - TreatmentFile (which should match the name where
%                           'Deisgn is saved')
%                       - DesignNumber
%                       - PlateShaking
%
%

document = com.mathworks.xml.XMLUtils.createDocument('Protocol');
protocol = document.getDocumentElement;

% Insert a comment to help document how the file was created.
protocol.appendChild(document.createTextNode(sprintf('\n   ')));
protocol.appendChild(document.createComment('Created by hpdd_exporter.m'));

% Constant for the non-ascii character mu.
MICRO = char(181);

% Add protocol-global settings.
create_text_children(protocol, ...
    {
    'Version'                   int2str(2);
    'VolumeUnit'                'nL';
    'ConcentrationUnit'         [MICRO 'M'];
    'MolarityConcentrationUnit' [MICRO 'M'];
    'MassConcentrationUnit'     'ng_mL';
    'ShakePerFluid'             logical2str(false);
    'ShakePlateDuration'        int2str(5);
    'ShakePerWell'              logical2str(true);
    'ShakeThresholdVolume'      int2str(100);
    'BackfillOrder'             'LastPriority';
    'BackfillNoDispense'        logical2str(0);
    } ...
    );

% Add Fluids list.
fluids = document.createElement('Fluids');
protocol.appendChild(fluids);
fluid_data = get_DesignDrugs(Designs);
% Map from drug name to string id for XML output.
fluid_ids = containers.Map;
for fluid_num = 1:length(fluid_data)
    fluid = document.createElement('Fluid');
    fluids.appendChild(fluid);
    % Use 0-based numbering for fluid IDs.
    id = int2str(fluid_num - 1);
    % Add to name->id map for later.
    fluid_ids(fluid_data(fluid_num).name) = id;
    fluid.setAttribute('ID', id);
    % The utility of RelatedID is not clear to me, but it was always -1 in
    % one sample file, and not present in another. -JLM
    %     related_id = int2str(-1);
    %     fluid.setAttribute('RelatedID', related_id);  %%%% discarded -MH
    % Concentration should be micromolar.
    conc_micromolar = fluid_data(fluid_num).stock_conc;
    ClassID = getVehicle(fluid_data(fluid_num).Vehicle);
    % Create fluid properties.
    create_text_children(fluid, ...
        {
        'Name'              fluid_data(fluid_num).name;
        'Concentration'     double2str(conc_micromolar);
        % It's not clear whether this is meant as the unit to use for
        % display within the D300 software, or something else, so we'll
        % just play it safe and use the same units we used up top.
        'ConcentrationUnit' [MICRO 'M'];
        % vehicule for the agent
        'ClassID' num2str(ClassID);
        } ...
        );
end

% Perform some checks on the Design data structures.
for design_num = 1:length(Designs)
    % Verify that the plate design has a standard format
    if mod(log2(Designs(design_num).plate_dims(1)),1) ~= 0 || ...
            (Designs(design_num).plate_dims(2)/Designs(design_num).plate_dims(1))~=1.5
        me = MException('ExportProtocol_D300:plate_layout_unknown', ...
            'Design %d is not a standard plate size', design_num);
        throw(me);
    end
    
    % Verify that all drug layout have the standard format
    for iD = 1:length(Designs(design_num).Drugs)
        if ~all(size(Designs(design_num).Drugs(iD).layout) == Designs(design_num).plate_dims)
            me = MException('ExportProtocol_D300:drug_layout_mismatch', ...
                'Drug %i for Design %d is not matching plate layout', iD, design_num);
            throw(me);
        end
    end
    
    % Verify that backfill maps are 2n x 3n (24/96/384/...) wells if it
    % exists. If we have to support different plate types in the future
    % this would need to be changed.
    if isfield(Designs(design_num), 'treated_wells') && ...
            ~all(size(Designs(design_num).treated_wells) == Designs(design_num).plate_dims)
        me = MException('ExportProtocol_D300:backfill_size_mismatch', ...
            'treated_wells for Design %d is not matching plate layout', design_num);
        throw(me);
    end
    
    if isfield(Designs(design_num), 'Vehicle') && ...
            ~all(size(Designs(design_num).Vehicle) == Designs(design_num).plate_dims)
        me = MException('ExportProtocol_D300:Vehicle_size_mismatch', ...
            'Vehicle for Design %d is not matching plate layout', design_num);
        throw(me);
    end
end

% load the table is a file name was passed
if ischar(t_plateinfo)
    t_plateinfo = tsv2table(t_plateinfo);
else
    assert(istable(t_plateinfo), ['Barcodes should be a table or a tsv file' ...
        ' with columns Barcode, DesignNumber'])
end

% Create Plates container.
plates = document.createElement('Plates');
protocol.appendChild(plates);
% NOTE What are the other Modes and do we need to support them?
create_text_children(plates, {'Mode' 'Concentration'});

% Check the different type of backfills
Vehicles = cell2mat(cellfun(@getVehicle, setdiff(unique([Designs.Vehicle]),''), ...
    'uniformoutput', false));

% Create Backfills container.
backfills = document.createElement('Backfills');
protocol.appendChild(backfills);
backfill = cell(1, length(Vehicles));
for iV = 1:length(Vehicles)
    backfill{iV} = document.createElement('Backfill');
    backfills.appendChild(backfill{iV});
    backfill{iV}.setAttribute('Type', 'ToMaxVolume');
    backfill{iV}.setAttribute('ClassID', num2str(Vehicles(iV)) );
    backfill_wells{iV} = document.createElement('Wells');
    backfill{iV}.appendChild(backfill_wells{iV});
    % FIXME For a six-plate experiment, saw two backfill elements each with half of
    %   the wells in it. Why?
end

assert(length(setdiff(unique(t_plateinfo.TreatmentFile),'-'))==1, ...
    ['Only one treatment file can be specified in the plate info file; ' ...
    'it should correspond to .mat file where the variable ''Designs'' is saved'])

t_trt_plates = t_plateinfo(~strcmp(t_plateinfo.TreatmentFile,'-'),:);

if ~isnumeric(t_trt_plates.DesignNumber)
    try
        t_trt_plates.DesignNumber = cellfun(@str2num,t_trt_plates.DesignNumber);
    catch
        disp('Plates with a TreatmentFile should have a numeric value for DesignNumber')
        disp(t_trt_plates(:,'DesignNumber'))
        error('Check plate info file')
    end
end

plate_cnt = 0;
for plate_num = 1:height(t_trt_plates)
    
    % Use barcode table's DesignNumber column as index into Design array.
    cur_design = Designs(t_trt_plates.DesignNumber(plate_num));
    if ~isfield(cur_design,'Drugs') || isempty(cur_design.Drugs) || ...
            all(reshape([cur_design.Drugs.layout],[],1)==0)
        warnprintf('Plate  %s  is ignored because no drug treatment was found', ...
            char(t_trt_plates.Barcode(plate_num)))
        continue
    end
    % if plates are skipped, the plate number in the hpdd file is not the
    % same as plate_num
    plate_cnt = plate_cnt + 1;
    
    plate_name = t_trt_plates.Barcode(plate_num);
    if isvariable(t_trt_plates, 'PlateShaking')
        PlateShaking = t_trt_plates.PlateShaking(plate_num);
    else
        PlateShaking = true;
    end
    % Convert volume from microliters to nanoliters.
    volume_nanoliters = cur_design.well_volume * 1e3;
    
    % creat the plate
    plate = document.createElement('Plate');
    plates.appendChild(plate);
    
    create_text_children(plate, ...
        {
        'PlateType'   sprintf('Default%i', numel(cur_design.Drugs(1).layout));
        'Rows'        int2str(size(cur_design.Drugs(1).layout,1));
        'Cols'        int2str(size(cur_design.Drugs(1).layout,2));
        'Name'        plate_name;
        'AssayVolume' double2str(volume_nanoliters);
        'DMSOLimit'   double2str(0.02);
        'AqueousLimit' double2str(0.05);
        'DontShake'   logical2str(~PlateShaking);
        % FIXME: AqueousLimit (may be an issue for old version of the hpdd
        % software (unchecked)
        } ...
        );
    wells = document.createElement('Wells');
    plate.appendChild(wells);
    for row = 1:size(cur_design.Drugs(1).layout,1)
        for column = 1:size(cur_design.Drugs(1).layout,2)
            % Create main well/fluid elements for drug treatments.
            well = document.createElement('Well');
            well.setAttribute('Row', int2str(row - 1));
            well.setAttribute('Col', int2str(column - 1));
            drugs = cur_design.Drugs;
            for drug_num = 1:length(drugs)
                % Values are already in the right units, micromolar.
                conc = drugs(drug_num).layout(row, column);
                if conc > 0
                    id = fluid_ids(drug2displayname(drugs(drug_num)));
                    fluid = document.createElement('Fluid');
                    fluid.setAttribute('ID', id);
                    fluid.setTextContent(double2str(conc));
                    well.appendChild(fluid);
                end
            end
            % Only add well element to document if not empty.
            if well.hasChildNodes
                wells.appendChild(well);
            end
            
            % Create backfill well element. If treated_wells is not present we
            % apply backfill to all wells. If it is present, we look up the
            % current well address in it to determine whether to apply backfill.
            if ~isfield(cur_design, 'treated_wells') || ...
                    cur_design.treated_wells(row,column)
                if isfield(cur_design, 'Vehicle') 
                    if ~isempty(cur_design.Vehicle{row,column})
                        iV = find(Vehicles == getVehicle(cur_design.Vehicle{row,column}));
                    else
                        continue
                    end
                else
                    iV = find(Vehicles == 0); % default is DMSO
                end
                backfill_well = document.createElement('Well');
                backfill_wells{iV}.appendChild(backfill_well);
                backfill_well.setAttribute('P', int2str(plate_cnt - 1));
                backfill_well.setAttribute('R', int2str(row - 1));
                backfill_well.setAttribute('C', int2str(column - 1));
            end
        end
    end
    % Even without randomization, the sample file had an empty Randomize tag.
    plate.appendChild(document.createElement('Randomize'));
end

xmlwrite([hpdd_filename '.hpdd'], document);

Write_DesignTreatment_summary([hpdd_filename '_summary.tsv'], Designs, t_trt_plates)

end


function create_text_children( element, data )
% Create a list of child nodes under element, using name/text pairs from
% the cell array in data.
document = element.getOwnerDocument;
for i = 1:size(data, 1)
    child = document.createElement(data{i,1});
    element.appendChild(child);
    child.setTextContent(data{i,2});
end
end




function name = drug2displayname(drug)
% Return the 'DrugName' and 'HMSLid' fields of drug, joined with a space.
name = drug.DrugName;
if length(drug.HMSLid) > 0
    name = [name ' ' drug.HMSLid];
end
end


function T = double2str(X)
% Format a double-precision float using 9 digits of precision, the same as
% the official D300 software appears to use.
T = sprintf('%.9g', X);
end


function T = logical2str(X)
% Format a logical as 'True' or 'False'.
if logical(X)
    T = 'True';
else
    T = 'False';
end
end

function ClassID = getVehicle(Vehicle)
if isnumeric(Vehicle)
    if ~ismember(Vehicle, 0:6)
        error('Unkonwn value (%f) for Vehicle', Vehicle)
    else
        ClassID = Vehicle;
    end
elseif ischar(Vehicle)
    switch lower(Vehicle)
        case {'dmso'}
            ClassID = 0;
        case {'aqueous + brij35' 'brij35'}
            ClassID = 2;2;
        case {'aqueous + brij35 + glycerol' 'brij35 + glycerol' 'brij35+glycerol'}
            ClassID = 5;
        case {'aqueous + tryton x100' 'tryton x100' 'tryton' 'trytonx100'}
            ClassID = 3;
        case {'aqueous + tryton x100 + glycerol' 'tryton x100 + glycerol' ...
                'tryton + glycerol' 'tryton+glycerol' 'trytonx100+glycerol'}
            ClassID = 6;
        case {'aqueous + tween 20' 'tween 20' 'tween' 'tween20'}
            ClassID = 1;
        case {'aqueous + tween 20 + glycerol' 'tween 20 + glycerol' ...
                'tween + glycerol' 'tween+glycerol' 'tween20+glycerol'}
            ClassID = 4;
        otherwise
            error('Unkonwn value (%s) for Vehicle', Vehicle)
    end
else
    error('Unkonwn class for Vehicle: %s', disp(Vehicle))
end
end
