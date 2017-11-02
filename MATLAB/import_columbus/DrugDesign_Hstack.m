function Design = DrugDesign_Hstack(Des1, Des2)
% Design = DrugDesign_Hstack(Des1, Des2)

%%
assert(length(Des1) == length(Des2))

Design = Des1;

for i=1:length(Des1)
    assert(Des1(i).well_volume == Des2(i).well_volume)
    Drugs = unique([{Des1(i).Drugs.DrugName} {Des2(i).Drugs.DrugName}])';
    DrugStruct = [ToRow(Des1(i).Drugs) ToRow(Des2(i).Drugs)];
    Design(i).treated_wells = [Des1(i).treated_wells Des2(i).treated_wells];
    Design(i).Vehicle = [Des1(i).Vehicle Des2(i).Vehicle];
    Design(i).plate_dims = size(Design(i).treated_wells);

    HMSLid = cell(length(Drugs),1);
    layout = cell(length(Drugs),1);
    stock_conc = cell(length(Drugs),1);
    Vehicle = cell(length(Drugs),1);
    for iD = 1:length(Drugs)
        HMSLid(iD) = unique({DrugStruct(strcmp({DrugStruct.DrugName}, Drugs{iD})).HMSLid});
        stock_conc{iD} = min([DrugStruct(strcmp({DrugStruct.DrugName}, Drugs{iD})).stock_conc]);
        assert(length(unique({DrugStruct(strcmp({DrugStruct.DrugName}, Drugs{iD})).Vehicle}))==1,...
            'inconsistent vehicles')
        Vehicle(iD) = unique({DrugStruct(strcmp({DrugStruct.DrugName}, Drugs{iD})).Vehicle});
        if ismember(Drugs{iD}, {Des1(i).Drugs.DrugName})
            layout1 = Des1(i).Drugs(strcmp(Drugs(iD), {Des1(i).Drugs.DrugName})).layout;
        else
            layout1 = zeros(Des1(i).plate_dims);
        end
        if ismember(Drugs{iD}, {Des2(i).Drugs.DrugName})
            layout2 = Des2(i).Drugs(strcmp(Drugs(iD), {Des2(i).Drugs.DrugName})).layout;
        else
            layout2 = zeros(Des2(i).plate_dims);
        end
        layout{iD} = [layout1 layout2];
    end

    Design(i).Drugs = struct('DrugName', Drugs, 'HMSLid', HMSLid, ...
        'stock_conc', stock_conc, 'layout', layout, 'Vehicle', Vehicle);


    Perturbations = unique([{Des1(i).Perturbations.Name} {Des2(i).Perturbations.Name}])';
    layout = cell(1,length(Perturbations));
     for iD = 1:length(Perturbations)
        if ismember(Perturbations{iD}, {Des1(i).Perturbations.Name})
            layout1 = Des1(i).Perturbations(strcmp(Perturbations(iD), {Des1(i).Perturbations.Name})).layout;
        else
            if iscell(Des2(i).Perturbations(strcmp(Perturbations(iD), ...
                    {Des2(i).Perturbations.Name})).layout)
                layout1 = cell(Des1(i).plate_dims);
            else
                layout1 = zeros(Des1(i).plate_dims);
            end
        end
        if ismember(Perturbations{iD}, {Des2(i).Perturbations.Name})
            layout2 = Des2(i).Perturbations(strcmp(Perturbations(iD), {Des2(i).Perturbations.Name})).layout;
        else
            if iscell(Des1(i).Perturbations(strcmp(Perturbations(iD), ...
                    {Des1(i).Perturbations.Name})).layout)
                layout2 = cell(Des2(i).plate_dims);
            else
                layout2 = zeros(Des2(i).plate_dims);
            end
        end
        layout{iD} = [layout1 layout2];
    end

    Design(i).Perturbations = struct('Name', ToRow(Perturbations), 'layout', ToRow(layout));


end
