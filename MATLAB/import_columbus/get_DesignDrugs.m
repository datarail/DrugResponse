
function drugs = get_DesignDrugs( design )
% Get list of unique drugs from design. Checks that stock_conc
% is the same for all same-named drugs.

drugs = struct('DrugName', {}, 'name', {}, 'HMSLid', {}, 'stock_conc', {}, 'Vehicle', {});
for design_num = 1:length(design)
    for drug_num = 1:length(design(design_num).Drugs)
        drug = design(design_num).Drugs(drug_num);

        DrugName = drug.DrugName;
        HMSLid = drug.HMSLid;
        if ~isempty(HMSLid)
            name = [DrugName ' ' HMSLid];
        else
            name = DrugName;
        end

        stock_conc = drug.stock_conc;
        if isfield(drug, 'Vehicle')
            Vehicle = drug.Vehicle;
        else
            Vehicle = 'DMSO'; % default is DMSO
        end
        match = find(strcmp(name, {drugs.name}));
        if ~isempty(match)
            if drugs(match).stock_conc ~= stock_conc
                me = MException(...
                    'hpdd_exporter:concentration_mismatch', ...
                    'Drug %s has different stock_conc values (%.1e %.1e)', ...
                    name, drugs(match).stock_conc, stock_conc ...
                );
                throw(me);
            end
            
            if drugs(match).Vehicle ~= Vehicle
                me = MException(...
                    'hpdd_exporter:concentration_mismatch', ...
                    'Drug %s has different Vehicles (%s - %s)', ...
                    name, drugs(match).Vehicle, Vehicle ...
                );
                throw(me);
            end
        else
            drugs(end+1) = struct( ...
                'DrugName', DrugName, ...
                'name', name, ...
                'HMSLid', HMSLid, ...
                'stock_conc', stock_conc, ...
                'Vehicle', Vehicle ...
            );
        end
    end
end
end
