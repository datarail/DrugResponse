function Design = DrugDesign_PutOnPlate(smallDesign, plate_dims, topcorner)
% Design = DrugDesign_PutOnPlate(smallDesign, plate_dims, topcorner)

%%


Design = smallDesign;



for i=1:length(Design)

    old_size = smallDesign(i).plate_dims;
    pos = false(plate_dims);
    pos(topcorner(1)+(0:(old_size(1)-1)), topcorner(2)+(0:(old_size(2)-1))) = true;

    Design(i).treated_wells = pos;
    Design(i).plate_dims = plate_dims;
    Design(i).Vehicle = repmat({''}, plate_dims);
    Design(i).Vehicle(topcorner(1)+(0:(old_size(1)-1)), topcorner(2)+(0:(old_size(2)-1))) = ...
        smallDesign(i).Vehicle;

    for iD = 1:length(Design(i).Drugs)
        Design(i).Drugs(iD).layout = zeros(plate_dims);
        Design(i).Drugs(iD).layout(pos(:)) = smallDesign(i).Drugs(iD).layout(:);
    end

    for iD = 1:length(Design(i).Perturbations)

        if iscell(smallDesign(i).Perturbations(iD).layout)
            Design(i).Perturbations(iD).layout = cell(plate_dims);
            Design(i).Perturbations(iD).layout(pos(:)) = smallDesign(i).Perturbations(iD).layout(:);
        else
            Design(i).Perturbations(iD).layout = zeros(plate_dims);
            Design(i).Perturbations(iD).layout(pos(:)) = smallDesign(i).Perturbations(iD).layout(:);
        end

    end
end
