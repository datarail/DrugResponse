function Write_DesignTreatment_summary(filename, Designs, t_plates)
% Write_D300_summary(filename, Designs, t_plates)

% load the table is a file name was passed
if ischar(t_plateinfo)
    t_plateinfo = tsv2table(t_plateinfo);
else
    assert(istable(t_plateinfo), ['Barcodes should be a table or a tsv file' ...
        ' with columns Barcode, DesignNumber'])
end

%%

Drugs = get_DesignDrugs( Designs(setdiff(unique(t_plates.DesignNumber),0)) );

output = ['Drug Name' 'HMSL id' 'Stock (mM)' ...
    strcat({t_plates.Barcode}, ' (ul)') 'D300 Total volume (ul)'];

Volume = zeros(length(Drugs), height(t_plates));

for iPlte = 1:height(t_plates)
    Didx = t_plates.DesignNumber(iPlte);
    if Didx==0
        continue
    end

    WellVolume = Designs(Didx).well_volume;
    for iDr=1:length(Designs(Didx).Drugs)

        Drug = Designs(Didx).Drugs(iDr).DrugName;
        Didx = find(strcmp({Drugs.name}, Drug));

        stock_conc = Drugs(Didx).stock_conc;

        Volume(Didx,iPlte) = sum(Designs(Didx).Drugs(iDr).layout(:)) ...
            *WellVolume/stock_conc;
    end
end

output(2:end,1) = {Drugs.name}';
output(2:end,2) = {Drugs.HMSLid}';
output(2:end,3) = {Drugs.stock_conc}';

for iDr=1:length(Drugs)
    output{iDr+1,end} = 1+ceil(sum(Volume(iDr,:)/100))/10;
end


%%
tsvwrite(filename, output, 'summary_drugs');


end
