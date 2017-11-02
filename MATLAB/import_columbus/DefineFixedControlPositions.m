function [ctrlidx, treated_wells] = DefineFixedControlPositions(plate_dims, ctrl_cnt)

fixedctrlpos = false(plate_dims);
treated_wells = true(plate_dims);

% specific position of the controls on the edges

assert(ctrl_cnt>=12, 'For controls on the edge, need at least 12 controls')

if ctrl_cnt<14 % only 6 controls on the edge
    disp('Only 6 controls on the edge, would be better with at least 14 controls in total')
    fixedctrlpos(round(1+((plate_dims(1)-1)*(0:3))/3), round(1+((plate_dims(2)-1)*(1:2))/3)) = true;
    fixedctrlpos(round(1+(plate_dims(1)-1)/2),round(1+((plate_dims(2)-1)*(0:3))/3)) = true;
else % 2 controls on each edge, 6 regularly spread in the middle.
    fixedctrlpos(round(1+((plate_dims(1)-1)*(0:3))/3), round(1+((plate_dims(2)-1)*(1:2))/3)) = true;
    fixedctrlpos(round(1+((plate_dims(1)-1)*(1:2))/3), round(1+((plate_dims(2)-1)*(0:3))/3)) = true;
end

if ctrl_cnt >= (sum(4*plate_dims)-16 +6) % remove two edges --> not treated
    fixedctrlpos([1 2 end-1 end], :) = true;
    fixedctrlpos(:, [1 2 end-1 end]) = true;
    treated_wells([1 2 end-1 end], :) = false;
    treated_wells(:, [1 2 end-1 end]) = false;

elseif ctrl_cnt >= (sum(2*plate_dims)-4 +6) % remove all edges --> not treated
    fixedctrlpos([1 end], :) = true;
    fixedctrlpos(:, [1 end]) = true;
    treated_wells([1 end], :) = false;
    treated_wells(:, [1 end]) = false;

elseif ctrl_cnt>=20 % put the corners as control (should be then discarded)
    fixedctrlpos([1 end], [1 end]) = true;
end

assert(all(size(fixedctrlpos)==plate_dims))
ctrlidx = find(fixedctrlpos);
