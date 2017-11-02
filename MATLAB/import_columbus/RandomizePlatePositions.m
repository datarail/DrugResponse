function DrugsStruct = RandomizePlatePositions(DrugsStruct, allTreatments, ...
    nWells, plate_dims, ctrlidx, ctrl_cnt, Seed)
% DrugsStruct = RandomizePlatePositions(DrugsStruct, allTreatments, ...
%     nWells, plate_dims, ctrlidx, ctrl_cnt, Seed+iR-1)

s = RandStream('mt19937ar','Seed',Seed);
RandStream.setGlobalStream(s);
if Seed>0
    idx = randperm(nWells); % randomize the order
else
    idx = 1:nWells;
end
idx(ctrlidx) = nWells+1; % put the 'fixed control' as control
order = sortidx(idx,'ascend'); % find the order for the treatment on the plate

for iD = 1:length(DrugsStruct)
    DrugsStruct(iD).layout(order(1:size(allTreatments,2))) = allTreatments(iD,:);
end

allDrugs = reshape([DrugsStruct.layout], [plate_dims length(DrugsStruct)]);
nDrugs = sum(allDrugs>0,3);

Treated = any(allDrugs>0,3);
assert(~any(Treated(ctrlidx)))
assert(sum(~Treated(:))==ctrl_cnt)
