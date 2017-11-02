function export_DrugDesign(Design1, filename, sheet)
% export_DrugDesign(Design1, filename, sheet)
%
%   export a Design structure (not an array) with standard fields for
%   treatment into tsv file with formated layout.
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO/drug/perturbation)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)

m = Design1.plate_dims(1);
n = Design1.plate_dims(2);

maxLabelName = 10; % length of drug name for layout sheet

total_output = cell(m+2,n+1);
letter_label = 'A':char('A'+m-1);
num_label = 1:n;

total_output(1+(1:m),1) = num2cell(letter_label);
total_output(1,2:end) = num2cell(num_label');

for i1 = 1:m
    for i2 = 1:n
        if ~Design1.treated_wells(i1,i2)
            total_output{i1+1,i2+1} = sprintf('Untreated well');
        end
    end
end

for i=1:length(Design1.Drugs)
    for i1 = 1:m
        for i2 = 1:n
            if Design1.Drugs(i).layout(i1,i2)>0
                if isempty(total_output{i1+1,i2+1})
                    total_output{i1+1,i2+1} = sprintf('%s, %.2fuM', ...
                        Design1.Drugs(i).DrugName(1:min(end,maxLabelName)), ...
                        Design1.Drugs(i).layout(i1,i2));
                else
                    total_output{i1+1,i2+1} = sprintf('%s\n%s, %.2fuM', ...
                        total_output{i1+1,i2+1}, Design1.Drugs(i).DrugName(1:min(end,maxLabelName)),...
                        Design1.Drugs(i).layout(i1,i2));
                end
            end
        end
    end

%     temp = {drugs_struct(i).name 'Stock (mM)=' drugs_struct(i).nominal_conc ...
%         'Volume (nl)=' drugs_struct(i).volume 'Well volume (ul)=' drugs_struct(i).well_volume*1e6 };
%     if isfield(drugs_struct(i),'Doses')
%         temp(2,1:(length(drugs_struct(i).Doses)+1)) = [{'Combo Doses (uM)='} ...
%             num2cell(drugs_struct(i).Doses)];
%     end
%     if isfield(drugs_struct(i),'SingleDoses')
%         temp(3,1:(length(drugs_struct(i).SingleDoses)+1)) = [{'Single Doses (uM)='} ...
%             num2cell(drugs_struct(i).SingleDoses)];
%     end
%     total_output(m+2+(i-1)*4+(1:3),1:size(temp,2)) = temp;

end


for i=1:length(Design1.Perturbations)

    for i1 = 1:m
        for i2 = 1:n
            total_output{i1+1,i2+1} = sprintf('%s\n%s=%s', ...
                total_output{i1+1,i2+1},Design1.Perturbations(i).Name(1:min(end,maxLabelName)), ...
                AnyToString(Design1.Perturbations(i).layout(i1,i2)));
        end
    end
end


if ~exist('filename','var')
    filename = 'temp.tsv';
end
if ~exist('sheet','var')
    sheet = 'default';
end

tsvwrite(filename, total_output, sheet)
