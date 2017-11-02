function t_data = PlateFilterByFocus(t_data, varargin)
% t_data = PlateFilterByFocus(t_data)



p = inputParser;
addParameter(p, 'pval', .01, @isnumeric)
addParameter(p, 'plotting', false, @islogical)

parse(p,varargin{:})
pval = p.Results.pval;
plotting = p.Results.plotting;
clear p


t_data.filtered = false(height(t_data),1);
t_data.deltafocus = NaN(height(t_data),1);

if ~all(isvariable(t_data, {'Column' 'Row'}))
    [Row,Column] = ConvertWellsToRowCol(cellstr(t_data.Well));
    t_data = [t_data table(Row,Column)];
end

plates = unique(t_data.Barcode);

for iP = 1:length(plates)
    pidx = t_data.Barcode == plates(iP);
    subt = t_data(pidx,:);

    [data, labels] = table_to_ndarray(subt, 'keyvars', {'Row' 'Column' 'Time'}, ...
        'outer', 1, 'valvars', 'focus');

    %%
    if plotting, get_newfigure(998,[50 550 900 300]), end

    for iC=1:length(labels{2}.Column)
        fprintf('.');
        for iR=1:length(labels{1}.Row)

            dist = (abs(labels{1}.Row-labels{1}.Row(iR)))*ones(1,height(labels{2})) + ...
                ones(height(labels{1}),1)*(.5+abs(labels{2}.Column-labels{2}.Column(iC)))';
            [test_R, test_C] = find(dist<=quantile(dist(:), 10/numel(dist)));
            focus = data(test_R, test_C, :);
            Wfocus = squeeze(data(iR,iC,:));

            preNdist = fitdist(focus(:), 'normal');
            Ndist = fitdist(focus( abs(preNdist.cdf(focus(:))-.5)<(.5-pval*2.5)), 'normal');
            Ndist = fitdist(focus( abs(Ndist.cdf(focus(:))-.5)<(.5-pval*2.5)), 'normal');

            if plotting,
                x = min(focus(:)):max(focus(:));
                prep = preNdist.pdf(x);
                p = Ndist.pdf(x);
                n = ksdensity(focus(:), x, 'width', diff(quantile(focus(:),[.48 .52])));
                nW = ksdensity(Wfocus, x, 'width', diff(quantile(focus(:),[.48 .52])));

                clf
                hold on
                plot(x, n/max(n), '-b')
                plot(x, nW/max(nW), '-b','linewidth',2)
                plot(x, prep/max(prep), '-r')
                plot(x, p/max(p), '-k','linewidth',2)
                lowerbound = Ndist.icdf(pval)-10;
                upperbound = Ndist.icdf(1-pval)+10;
                plot([lowerbound*[1 1] upperbound*[1 1]], [.5 0 0 .5], '-k')
                if any( (Wfocus<lowerbound) | (Wfocus>upperbound))
                    plot(squeeze(data(iR,iC,(Wfocus<lowerbound) | (Wfocus>upperbound))), .5, '*')
                    pause
                end
            end
            %%

            idx = find(pidx & t_data.Row==labels{1}.Row(iR) & ...
                t_data.Column==labels{2}.Column(iC));
            [~, idxT] = ismember(labels{3}.Time, t_data.Time(idx));


            t_data.filtered(idx(idxT)) = (Wfocus<lowerbound) | (Wfocus>upperbound);
            t_data.deltafocus(idx(idxT)) = Wfocus-Ndist.mu;
        end
        %         if any(abs(delta(:))>1)
        %             [filt_Row, filt_Col] = find(abs(delta)>1); true;
        %         end

    end
    fprintf('\n');
end
