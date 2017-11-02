function [Design, t_barcode] = hpdd_importer(hpdd_filename)
% [Design, t_barcode] = hpdd_importer(hpdd_filename)
%   Read a hpdd file and convert the information in an array of Design structure
%   with standard fields for treatment. Convert the metadata in a table with
%   barcodes and corresponding treatment design.
%
%   hpdd_filename : path and file name to a treatment saved by the D300
%   Design :        array of design structures with the following fields:
%                       - plate_dims (plate dimension)
%                       - treated_wells (wells treated with DMSO)
%                       - well_volume (in uL)
%                       - Drugs (structure with DrugName, HMSL_id, stock_conc
%                           and layout - concentration given in uM)
%   t_barcode :     table with columns:
%                       - Barcode
%                       - TreatmentFile (input)
%                       - DesignNumber
%

%% get the data

protocol = XMLElement.parse(hpdd_filename);

%% get the main fields of the xml structure and perform a few checks
DMSObackfill = protocol.find('Backfills').iter('Backfill');
plates = protocol.find('Plates').iter('Plate');
Drugstruct = protocol.find('Fluids').iter('Fluid');

% unit consistency
assert(strcmp(protocol.find('ConcentrationUnit').text, ...
    protocol.find('MolarityConcentrationUnit').text), ...
    'Mismatch between the units for ''ConcentrationUnit'' and ''MolarityConcentrationUnit''');


%% properties of the drugs

DrugNames = cell(length(Drugstruct),1);
stock_conc = cell(length(Drugstruct),1);
DrugIdxes = NaN(length(Drugstruct),1);

for i = 1:length(Drugstruct)

    % check for the proper ordering (necessary for indexing)
    DrugIdxes(i) = str2double(Drugstruct(i).attributes(strcmp({Drugstruct(i).attributes.name},'ID')).value);

    DrugNames{i} = Drugstruct(i).find('Name').text;

    % convert to uM for the stock concentration
    if Drugstruct(i).find('ConcentrationUnit').text(1) == 'n' % nM
        Concentration_conversion = 1e-3;
    elseif Drugstruct(i).find('ConcentrationUnit').text(1) == 181 % uM
        Concentration_conversion = 1;
    elseif Drugstruct(i).find('ConcentrationUnit').text(1) == 'm' % mM
        Concentration_conversion = 1e3;
    else
        error('Issue with the unit of the stock drug concentration: %s', ...
            Drugstruct(i).find('ConcentrationUnit').text)
    end

    stock_conc{i} = str2double(Drugstruct(i).find('Concentration').text) ...
        *Concentration_conversion;

end

[DrugNames, HMSLids] = splitHMSLid(DrugNames);

Drugs = struct('DrugName',DrugNames, 'stock_conc',stock_conc, ...
    'layout', []);

if any(~cellfun(@isempty, HMSLids))
    for i=1:length(Drugs)
        Drugs(i).HMSLid = HMSLids{i};
    end
end

%% get the general properties of the plates

plate_dims = cell(length(plates),1);
treated_wells = plate_dims;
well_volume = plate_dims;   % will be transformed in uL
Barcode = plate_dims;
treated_wells = plate_dims;


if protocol.find('VolumeUnit').text(1) == 'n' % nL
    volume_conversion = 1e-3;
elseif protocol.find('VolumeUnit').text(1) == 181 % uL
    volume_conversion = 1;
else
    error('Issue with the unit of the well volume: %s', ...
        protocol.find('VolumeUnit').text)
end


for iP = 1:length(plates)
    plate_dims{iP} = [str2double(plates(iP).find('Rows').text) ...
        str2double(plates(iP).find('Cols').text)];
    well_volume{iP} = str2double(plates(iP).find('AssayVolume').text) ...
        *volume_conversion;
    Barcode{iP} = plates(iP).find('Name').text;

    treated_wells{iP} = false(plate_dims{iP});
end

for iB = 1:length(DMSObackfill)
    Wells = DMSObackfill(iB).find('Wells').iter('Well');
    for iW = 1:length(Wells)
        well = Wells(iW);
        row = str2double(well.get('R')) + 1;
        col = str2double(well.get('C')) + 1;
        iP = str2double(well.get('P')) + 1;

        assert(iP<=length(plates))
        assert(row<=plate_dims{iP}(1))
        assert(col<=plate_dims{iP}(2))

        treated_wells{iP}(row, col) = true;
    end
end


Design = struct('plate_dims', plate_dims, 'treated_wells', treated_wells, ...
    'well_volume', well_volume, 'Drugs', Drugs);


%% assign the concentration for each drug

DMSOwarning = true;

for iP = 1:length(plates)

    Wells = plates(iP).find('Wells').iter('Well');
    % randomization
    randomized = ~isempty(plates(iP).find('Randomize')) & ...
        ~isempty(plates(iP).find('Randomize').iter('Well'));
    if randomized
        rWells = plates(iP).find('Randomize').iter('Well');
        assert(length(rWells)== ...
            max(length(Wells),sum(treated_wells{iP}(:))))
        Randomization = NaN(length(rWells),4);
        for i=1:length(rWells)
            assert(all(strcmp({rWells(i).attributes.name}, {'C' 'R' 'ToC' 'ToR'})))
            Randomization(i,:) = cellfun(@str2double, {rWells(i).attributes.value});
        end
        Randomization = Randomization+1;
    end

    for i=1:length(Design(iP).Drugs)
        Design(iP).Drugs(i).layout = zeros(Design(iP).plate_dims);
    end
    usedDrug = false(length(Design(iP).Drugs), 1);

    for iW = 1:length(Wells)
        well = Wells(iW);
        row = str2double(well.get('Row')) + 1;
        col = str2double(well.get('Col')) + 1;

        if randomized
            idx = find(all(Randomization(:,3)==col & Randomization(:,4)==row,2));
            assert(length(idx)==1)
            col = Randomization(idx,1);
            row = Randomization(idx,2);
        end

        assert(row<=Design(iP).plate_dims(1))
        assert(col<=Design(iP).plate_dims(2))

        fluids = well.iter('Fluid');
        for iD = 1:length(fluids)
            Didx = find( str2double(fluids(iD).get('ID')) == DrugIdxes);
            assert(Didx>0 && Didx<=length(Design(iP).Drugs))
            usedDrug(Didx) = true;
            Design(iP).Drugs(Didx).layout(row, col) = str2double(fluids(iD).text);
            if ~(Design(iP).treated_wells(row, col))
                if DMSOwarning
                    warning('Some wells have no DMSO backfill')
                end
                Design(iP).treated_wells(row, col) = true;
            end
        end
    end

    % remove the drugs that are not used in that particular design
    Design(iP).Drugs = Design(iP).Drugs(usedDrug);
    assert(all(cellfun(@(x)any(any(x>0)), {Design(iP).Drugs.layout})))

end

%% find the replicates and match the bar codes.

DesignNumber = (1:length(Barcode))';
DesignToRemove = false(length(Barcode),1);


for iP = 2:length(plates)
    for iP2 = 1:(iP-1)
        if isequal(Design(iP2), Design(iP))
            DesignToRemove(iP) = true;
            DesignNumber(iP) = iP2;
            break
        end
    end
end

DesignIdx = DesignNumber(~DesignToRemove);
assert(all(ismember(DesignNumber, DesignIdx)))

Design = Design(DesignIdx);
DesignNumber = arrayfun(@(x) find(x==DesignIdx), DesignNumber);

%% table with barcodes
[~,fname,fext] = fileparts(hpdd_filename);
TreatmentFile = repmat({[fname fext]}, length(plates),1);

t_barcode = table(Barcode, TreatmentFile, DesignNumber);
