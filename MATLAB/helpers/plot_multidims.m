function [a, hl] = plot_multidims(t_data, varargin)

%  [a, hl] = plot_multidims(t_data, varargin)
%       inputs: xplotkey, yplotkey, xaxiskey, yaxiskey, colorkey,
%       xtransform, ytransform, axischanges, mean_SEM, xspacing, yspacing,
%       yval_lines, plotcolors
%

Generate_Plotting_parameters

noyplot = 'noyplot';
nocplot = 'nocplot';

p = inputParser;
addParameter(p, 'xplotkey', @(x) ischar(x) || iscellstr(x))
addParameter(p, 'yplotkey', noyplot, @(x) ischar(x) || iscellstr(x))
addParameter(p, 'xaxiskey', @isstr)
addParameter(p, 'yaxiskey', @isstr)
addParameter(p, 'colorkey', nocplot, @isstr)
addParameter(p, 'xtransform', @(x)x, @(x)isa(x,'function_handle'));
addParameter(p, 'ytransform', @(x)x, @(x) isa(x,'function_handle') || ...
    all(cellfun(@(y) isa(y,'function_handle'), x)))
addParameter(p, 'axischanges', @(x) set(x,'fontsize',6), @(x) isa(x,'function_handle'))
addParameter(p, 'mean_SEM', true, @islogical)
addParameter(p, 'plotcolors', jet(20), @(x) isnumeric(x) || isa(x,'function_handle'))
addParameter(p, 'xspacing', .03, @isnumeric)
addParameter(p, 'yspacing', .07, @isnumeric)
addParameter(p, 'yval_lines', [0 1], @isnumeric)

parse(p,varargin{:})
p = p.Results;

if ischar(p.xplotkey)
    p.xplotkey = {p.xplotkey};
end
if isa(p.plotcolors,'function_handle') && ~strcmp(p.colorkey, nocplot)
    p.plotcolors = p.plotcolors(length(unique(t_data.(p.colorkey))));
end

xplotkeys = unique(t_data(:,p.xplotkey));
if strcmp(p.yplotkey, noyplot)
    yplotkeys = table(0,'variablename',{'noyplot'});
    t_data = [t_data table(zeros(height(t_data),1), 'variablenames', {noyplot})];
elseif strcmp(p.yplotkey, 'xkeys')
    yplotkeys = table(1,'variablename',{'noyplot'});
    p.yplotkey = noyplot;
    t_data = [t_data table(ones(height(t_data),1), 'variablenames', {noyplot})];
else
    yplotkeys = unique(t_data(:,p.yplotkey));
end
if ischar(p.yplotkey)
    p.yplotkey = {p.yplotkey};
end
if strcmp(p.colorkey, nocplot)
    colorkeys = 1;
    t_data = [t_data table(ones(height(t_data),1), 'variablenames', {nocplot})];
else
    colorkeys = unique(t_data.(p.colorkey));
end
if ischar(p.yaxiskey)
    p.yaxiskey = {p.yaxiskey};
end
if isa(p.ytransform,'function_handle')
    p.ytransform = repmat({p.ytransform},length(p.yaxiskey),1);
end
%
t_data = TableToCategorical(t_data, [p.xplotkey, p.yplotkey, {p.colorkey}]);

%%

if ~isvariable(yplotkeys,noyplot) || yplotkeys.noyplot~=1
    nCols = size(xplotkeys,1);
    nRows = size(yplotkeys,1);
else
    nCols = ceil(sqrt(size(xplotkeys,1)));
    nRows = ceil(size(xplotkeys,1)/nCols);
end

xspacing = p.xspacing;
axis_width = (1-(nCols+1)*xspacing)/nCols;
assert(axis_width>0, 'Too many columns, reduce xspacing')
yspacing = p.yspacing;
axis_height = (1-(nRows+1)*yspacing)/nRows;
assert(axis_height>0, 'Too many rows, reduce yspacing')

linetypes = {'-' ':' '--'};
for iyp = 1:nRows
    for ixp=1:nCols


        if isvariable(yplotkeys,'noyplot')
            if yplotkeys.noyplot==0
                xidx = ixp;
            elseif yplotkeys.noyplot==1
                xidx = (iyp-1)*nCols + ixp;
                if xidx>size(xplotkeys,1)
                    break
                end
            end
            yidx = 1;
        else
            xidx = ixp;
            yidx = iyp;
        end

        a(iyp, ixp) = get_newaxes([xspacing*1.5+(ixp-1)*(axis_width+xspacing) ...
            yspacing*1.5+(nRows-iyp)*(axis_height+yspacing) axis_width axis_height],1,...
            'fontsize',6);


        if isvariable(yplotkeys,'noyplot')
            title([strjoin(p.xplotkey,',') '=' strjoin(table2cellstr(xplotkeys(xidx,:),0),',')])
        else
            title({[strjoin(p.xplotkey,',') '=' strjoin(table2cellstr(xplotkeys(xidx,:),0),',') ';'];
                [strjoin(p.yplotkey,',') '=' strjoin(table2cellstr(yplotkeys(yidx,:),0),',')]})
        end


        xvals = t_data.(p.xaxiskey)(eqtable(t_data(:,p.xplotkey),xplotkeys(xidx,:)) & ...
            eqtable(t_data(:,p.yplotkey),yplotkeys(yidx,:)));
        if isempty(xvals), continue, end
        xvals = p.xtransform(xvals);
        for i=1:length(p.yval_lines)
            plot([min(xvals(:)) max(xvals(:))], p.yval_lines(i)*[1 1], ...
                '-', 'color', [.8 .8 .8])
        end

        for iC = 1:length(colorkeys)
            subt = t_data(eqtable(t_data(:,p.xplotkey),xplotkeys(xidx,:)) & ...
                eqtable(t_data(:,p.yplotkey),yplotkeys(yidx,:)) & ...
                t_data.(p.colorkey)==colorkeys(iC),[{p.xaxiskey} p.yaxiskey]);
            if isempty(subt)
                continue
            end

            subt.(p.xaxiskey) = p.xtransform(subt.(p.xaxiskey));
                for iyak = 1:length(p.yaxiskey)
                    subt.(p.yaxiskey{iyak}) = p.ytransform{iyak}(subt.(p.yaxiskey{iyak}));
                end
            subt = sortrows(subt, p.xaxiskey);

            if ~p.mean_SEM
                for iyak = 1:length(p.yaxiskey)
                    temph = plot(subt.(p.xaxiskey), subt.(p.yaxiskey), ['.' linetypes{iyak}], ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    if iyak==1, h(iC) = temph;end
                end

            else
                y_mean = collapse(subt, @mean, 'keyvars', {p.xaxiskey}, 'valvars', p.yaxiskey);
                y_SEM = collapse(subt, @SEM,  'keyvars', {p.xaxiskey}, 'valvars', p.yaxiskey);
                for iyak = 1:length(p.yaxiskey)
                    if height(y_SEM) == height(subt)
                        temph = plot(subt.(p.xaxiskey), subt.(p.yaxiskey{iyak}), ['.' linetypes{iyak}], ...
                            'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    else
                        temph = errorbar(y_mean.(p.xaxiskey), y_mean.(p.yaxiskey{iyak}),...
                            y_SEM.(p.yaxiskey{iyak}), ['.' linetypes{iyak}], ...
                            'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                    end
                    if iyak==1, h(iC) = temph;end
                end
            end
        end

        p.axischanges(gca)
        if iyp==nRows, xlabel(gca,p.xaxiskey,'fontweight','bold'), end
        if ixp==1, ylabel(gca, strjoin(p.yaxiskey),'fontweight','bold'), end

        if ((ixp==nCols && iyp==nRows) || xidx==size(xplotkeys,1)) && ~strcmp(p.colorkey, nocplot)            
            hl = legend(h(ishandle(h)), strcat(p.colorkey, '=', AnyToString(colorkeys(ishandle(h)))), ...
                'fontsize',6, 'orientation', 'horizontal');
            set(hl, 'position', [.01 .005 .98 .03])
        end
    end
end

if nargout==0
 clear a
end
