
function t_processed_t0 = TimeCourseExtrapolate_t0(t_processed)
%
% t_mean_t0 = TimeCourseExtrapolate_t0(t_processed);
%
%

t_plate_t0 = unique(t_processed(:,'Barcode'));
t_plate_t0.Day0Cnt(:) = NaN;
for iP = 1:height(t_plate_t0)

    subt = sortrows(unique(t_processed(t_processed.Barcode==t_plate_t0.Barcode(iP),...
        {'Barcode' 'Ctrlcount' 'Time'})),3);

    cnt_t0 = exp(log(subt.Ctrlcount(1)) - subt.Time(1)*...
        (log(subt.Ctrlcount)-log(subt.Ctrlcount(1)))./(subt.Time-subt.Time(1)));
    t_plate_t0.Day0Cnt(iP) = mean(cnt_t0([2 3 3 4 4 5]));

end

%%
t_processed_t0 = leftjoin(t_processed(t_processed.DrugName~='-',:), t_plate_t0);

t_processed_t0 = [t_processed_t0 array2table([ ...
    (t_processed_t0.Cellcount-t_processed_t0.Day0Cnt)./(t_processed_t0.Ctrlcount-t_processed_t0.Day0Cnt) ...
    2.^(log2(t_processed_t0.Cellcount./t_processed_t0.Day0Cnt)./log2(t_processed_t0.Ctrlcount./t_processed_t0.Day0Cnt))-1], ...
    'variablenames', {'RelGrowth' 'nRelGrowth'})];
