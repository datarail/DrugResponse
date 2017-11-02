function [biased, edge_res, col_res, row_res] = TestPlateBias(t_data, BiasCutoff, plotting)
% [biased, edge_res, col_res, row_res] = TestPlateBias(t_data, BiasCutoff, plotting)
%       t_data:     need variables 'Well' or 'Row'/'Column'; all other
%                       columns will be checked against.
%
%       BiasCutoff: [minimal drop cutoff      p-value cutoff];
%                       default = [.1 .01]
%                   can be a 3x2 matrix for each condition:
%                           edge, column, row
%
%       biased:     [edge_res, col_res, row_res] < cutoff & significant
%
%       xxx_res:    [mean(vals) std(vals) p];
%
%

if ~exist('plotting','var')
    plotting = .5;
end


if ~exist('BiasCutoff','var') || isempty('BiasCutoff')
    BiasCutoff = [1;1;1]*[.1 .01];
elseif all(size(BiasCutoff)==[1 2])
    BiasCutoff = [1;1;1]*BiasCutoff;
elseif ~all(size(BiasCutoff)==[3 2])
    error('Wrong size for BiasCutoff input, needs to be 1x2 or 3x2')
end

if ~all(isvariable(t_data, {'Column' 'Row'}))
    [Row,Column] = ConvertWellsToRowCol(cellstr(t_data.Well));
    t_data = [t_data table(Row,Column)];
    t_data.Well = [];
end

if isvariable(t_data,'Barcode')
    barcode = unique(t_data.Barcode);
    assert(length(barcode)==1, 'More than one barcode in t_data')
    barcode = char(barcode);
else
    barcode = '';
end

if isvariable(t_data,'Time')
    Time = unique(t_data.Time);
    assert(length(Time)==1, 'More than one time point in t_data')
    barcode = [barcode ' - ' num2str(Time) 'h'];
end

assert(height(t_data)==height(unique(t_data(:,{'Row' 'Column'}))), ...
    'Multiple values for the same well in t_data');

%%
maxrow = 2^ceil(log2(max(t_data.Row)));
maxcol = 3*2^ceil(log2(max(t_data.Column)/3));
t_data.DistEdge = min([(t_data.Row) (maxrow-t_data.Row+1) ...
    (t_data.Column) (maxcol-t_data.Column+1)],[],2);

VarToTest = setdiff(varnames(t_data), {'Well' 'Row' 'Column' 'DistEdge' 'Barcode' 'Time'});

[plate, labels] = table2ndarray(t_data, 'keyvars', {'Row' 'Column'}, ...
    'outer', 1, 'valvars', VarToTest);


for iv = 1:length(VarToTest)

    meanvals = mean(t_data.(VarToTest{iv}));
    relvals = (t_data.(VarToTest{iv})-meanvals)/meanvals;

    %% test for biased based on the distance from the edge
    edge_res = NaN(min(6, max(t_data.DistEdge)), 3);
    for i=1:min(6, max(t_data.DistEdge))
        % takes two ranges at the same time
        vals = relvals(abs(t_data.DistEdge-i)<2);
        nulldist = relvals(t_data.DistEdge>(i+1));

        if length(vals)<30 || length(nulldist)<30
            continue
        end

        [~,p] = ttest2(vals, nulldist);

        edge_res(i,:) = [mean(vals) std(vals) p];
    end

    %% test for biased based on the column (set of 3 columns)

    col_res = NaN(maxcol,3);
    for i=1:maxcol
        vals = relvals( abs(t_data.Column-i)<2 );
        nulldist = relvals( abs(t_data.Column-i)>=2 );

        if length(vals)<30 || length(nulldist)<30, continue, end

        [~,p] = ttest2(vals, nulldist);

        col_res(i,:) = [mean(vals) std(vals) p];
    end

    %% test for biased based on the column (set of 2 columns)

    row_res = NaN(maxrow,3);
    for i=1:maxrow
        vals = relvals( abs(t_data.Row-i)<2 );
        nulldist = relvals( abs(t_data.Row-i)>=2 );

        if length(vals)<30 || length(nulldist)<30, continue, end

        [~,p] = ttest2(vals, nulldist);

        row_res(i,:) = [mean(vals) std(vals) p];
    end

    %% plot results
    biased = [any(edge_res(:,1)<-BiasCutoff(1,1) & edge_res(:,3)<BiasCutoff(1,2)) ...
        any(col_res(:,1)<-BiasCutoff(2,1) & col_res(:,3)<BiasCutoff(2,2)) ...
        any(row_res(:,1)<-BiasCutoff(3,1) & row_res(:,3)<BiasCutoff(3,2))];

    %
    if plotting>=1 || (plotting==.5 && any(biased))
        figure(999)
        clf
        set(gcf, 'filename', ['Test_bias' barcode '_' VarToTest{iv} '.pdf'])
        get_newaxes([.05 .1 .55 .85])
        imagesc(labels{2}.Column, labels{1}.Row, plate(:,:,iv), ...
            [0 max(max(plate(:,:,1)))])
        xlim([.5 maxcol+.5])
        ylim([.5 maxrow+.5])

        title([barcode ' - ' VarToTest{iv}], 'fontsize', 8, 'fontweight', 'bold', ...
            'interpreter', 'none')
        set(gca,'fontsize',6, 'ytick', labels{1}.Row, 'xtick', labels{2}.Column(1:2:end))


        get_newaxes([.65 .7 .31 .24],1)
        h = xyerrorbars(1:size(edge_res,1), [], edge_res(:,1), edge_res(:,2));
        set(h(1),'linewidth',1.5, 'color','k')
        plot([0 size(edge_res,1)+[1 1] 0], -BiasCutoff(1,1)*[1 1 0 0], 'k');
        for i=1:size(edge_res,1)
            if edge_res(i,3)<BiasCutoff(1,2)
                text(i, .1, '*', 'fontsize', 7, 'fontweight', 'bold', ...
                    'horizontalalign', 'center')
                if edge_res(i,1)<-BiasCutoff(1,1)
                    plot(i, edge_res(i,1), 'ok')
                end
            end
        end
        xlim([.5 size(edge_res,1)+.5])
        ylim([-.7 .5])
        title('Edge bias', 'fontsize', 8, 'fontweight', 'bold')
        set(gca,'fontsize',6)

        get_newaxes([.65 .37 .31 .24],1)
        h = xyerrorbars(1:size(col_res,1), [], col_res(:,1), col_res(:,2));
        set(h(1),'linewidth',1.5, 'color','k')
        plot([0 size(col_res,1)+[1 1] 0], -BiasCutoff(2,1)*[1 1 0 0], 'k');
        for i=1:size(col_res,1)
            if col_res(i,3)<BiasCutoff(2,2)
                text(i, .1, '*', 'fontsize', 7, 'fontweight', 'bold', ...
                    'horizontalalign', 'center')
                if col_res(i,1)<-BiasCutoff(2,1)
                    plot(i, col_res(i,1), 'ok')
                end
            end
        end
        xlim([.5 size(col_res,1)+.5])
        ylim([-.7 .5])
        title('Column bias', 'fontsize', 8, 'fontweight', 'bold')
        set(gca,'fontsize',6)

        get_newaxes([.65 .04 .31 .24],1)
        h = xyerrorbars(1:size(row_res,1), [], row_res(:,1), row_res(:,2));
        set(h(1),'linewidth',1.5, 'color','k')
        plot([0 size(row_res,1)+[1 1] 0], -BiasCutoff(3,1)*[1 1 0 0], 'k');
        for i=1:size(row_res,1)
            if row_res(i,3)<BiasCutoff(3,2)
                text(i, .1, '*', 'fontsize', 7, 'fontweight', 'bold', ...
                    'horizontalalign', 'center')
                if row_res(i,1)<-BiasCutoff(3,1)
                    plot(i, row_res(i,1), 'ok')
                end
            end
        end
        xlim([.5 size(row_res,1)+.5])
        ylim([-.7 .5])
        title('Row bias', 'fontsize', 8, 'fontweight', 'bold')
        set(gca,'fontsize',6)

        if plotting==1
            pause
        elseif plotting==2
            drawnow
        else
            pause(1)
        end

    end

end
