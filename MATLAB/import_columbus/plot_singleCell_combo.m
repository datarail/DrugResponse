function plot_singleCell_combo(OutputMx, OutputDist, Concs, Distbins, pos)

if ~exist('pos','var')
    pos = [.15 .15 .8 .8];
end


lC1 = length(Concs{1});
lC2 = length(Concs{2});
lDb = length(Distbins);

get_newaxes(pos,1)

imagesc(1:lC2, 1:lC1, OutputMx(:,:,1), [-.5 1.5])
cmap = [1 1 1;
    (.31:.01:.8)' (.01:.01:.5)'*[1 1]
    min((.81:.01:1.8),1)' ((.505:.005:1)'*[1 1]).^2
    ((1:-.01:.51)'*[1 1]).^2 ones(50,1)];
colormap(cmap)

xlim([.5 lC2+.5])
ylim([.5 lC1+.5])

for iC1 = 1:lC1
    for iC2 = 1:lC2
        xpos = iC2-.46;
        ypos = iC1-.4;
        if ~isempty(OutputDist{iC1, iC2})
            plot( xpos+[0 0 .92], ypos+[.85 0 0], '-k')
            plot( xpos+[1 1]*find(Distbins==0)/lDb, ypos+[-.05 .85], ':k')
            plot( xpos+.9*(1:lDb)/lDb, ypos+.8*OutputDist{iC1, iC2}./...
                (ones(lDb,1)*max(OutputDist{iC1, iC2})), '-', ...
                'color', [.4 .6 .4])
            plot( xpos+.9*(1:lDb)/lDb, ypos+.8*mean(OutputDist{iC1, iC2},2)./...
                (ones(lDb,1)*max(mean(OutputDist{iC1, iC2},2))), '-k')
            plot( [1 1]*xpos+.9*find(cumsum(mean(OutputDist{iC1, iC2},2))>...
                .5*sum(mean(OutputDist{iC1, iC2},2)),1,'first')/lDb, ...
                ypos+[-.05 .85], '-k')
        end
    end
end

set(gca, 'xtick', 1:lC2, 'xticklabel', num2cellstr(Concs{2},'%.2g'), 'ytick', 1:lC1, ...
    'yticklabel', num2cellstr(Concs{1},'%.2g'), 'fontsize', 6)
