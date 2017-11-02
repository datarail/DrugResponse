function [Designs, Designs_DMSO] = DrugDesign_AddDyeDMSO(Designs, name, conc, stock)
% [Designs, Designs_DMSO] = DrugDesign_AddDyeDMSO(Designs, name, conc, stock)

if ischar(name), name = {name}; end

assert(length(name)==length(conc))
assert(length(name)==length(stock))

Designs_DMSO = Designs;


for i=1:length(Designs)
    for iD = 1:length(name)

        Designs_DMSO(i).Drugs(end+1) = struct('DrugName', name{iD}, 'HMSLid', '-', ...
            'stock_conc', stock(iD), 'layout', conc(iD)*Designs(iD).treated_wells);

        Designs(i).Perturbations(end+1) = struct('Name', name{iD},  ...
            'layout', conc(iD)*Designs(iD).treated_wells);

    end
end
