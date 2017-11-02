
function Design = SetDrugVehicle(Design)
% Design = SetDrugVehicle(Design)

% assign the vehicle for the treated wells
[nV, Vehicles] = hist_str({Design.Drugs.Vehicle});
trt = Design.treated_wells;
Design.Vehicle = cell(size(trt));
if length(Vehicles)>1
    SelectVec = Vehicles{argmax(nV)};
    warnprintf('multiple vehicles: setting ALL controls to %s', SelectVec)
else
    SelectVec = Vehicles{1};
end
for iV = 1:length(Vehicles)
    trt_veh = reshape(...
        [Design.Drugs(strcmp(Vehicles{iV}, {Design.Drugs.Vehicle})).layout], ...
        Design.plate_dims(1), Design.plate_dims(2), []);
    idx = any(trt_veh>0,3);
    warnassert(all(cellfun(@isempty, Design.Vehicle(idx))), ...
        'Some wells have drugs with different vehicles');
    Design.Vehicle(idx) = Vehicles(iV);
end
assert(all(cellfun(@isempty, Design.Vehicle(~trt))), ...
    'Inconsistency between treatments and controls');
Design.Vehicle(cellfun(@isempty, Design.Vehicle)) = {SelectVec};
