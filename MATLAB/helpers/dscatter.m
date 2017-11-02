function [hAxes, maxF, edges1, edges2] = dscatter(X,Y, varargin)
% DSCATTER creates a scatter plot coloured by density.
%
%   DSCATTER(X,Y) creates a scatterplot of X and Y at the locations
%   specified by the vectors X and Y (which must be the same size), colored
%   by the density of the points.
%
%   DSCATTER(...,'MARKER',M) allows you to set the marker for the
%   scatter plot. Default is 's', square.
%
%   DSCATTER(...,'MSIZE',MS) allows you to set the marker size for the
%   scatter plot. Default is 10.
%
%   DSCATTER(...,'FILLED',false) sets the markers in the scatter plot to be
%   outline. The default is to use filled markers.
%
%   DSCATTER(...,'PLOTTYPE',TYPE) allows you to create other ways of
%   plotting the scatter data. Options are "surf','mesh' and 'contour'.
%   These create surf, mesh and contour plots colored by density of the
%   scatter data.
%
%   DSCATTER(...,'BINS',[NX,NY]) allows you to set the number of bins used
%   for the 2D histogram used to estimate the density. The default is to
%   use the number of unique values in X and Y up to a maximum of 200.
%
%   DSCATTER(...,'SMOOTHING',LAMBDA) allows you to set the smoothing factor
%   used by the density estimator. The default value is 20 which roughly
%   means that the smoothing is over 20 bins around a given point.
%
%   DSCATTER(...,'LOGX',true,'LOGY',true) log10 the x/y values.
%
%   DSCATTER(...,'LOGC',true) uses a log scale for the color.
%
%   DSCATTER(...,'maxF',maxF) uses maxF to scale colors to the top value (cap to 1).
%
%   Examples:
%
%       [data, params] = fcsread('SampleFACS');
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1)
%       % Add contours
%       hold on
%       dscatter(data(:,1),10.^(data(:,2)/256),'log',1,'plottype','contour')
%       hold off
%       xlabel(params(1).LongName); ylabel(params(2).LongName);
%       
%   See also FCSREAD, SCATTER.

% Copyright 2003-2004 The MathWorks, Inc.
% $Revision:  $   $Date:  $

% Reference:
% Paul H. C. Eilers and Jelle J. Goeman
% Enhancing scatterplots with smoothed densities
% Bioinformatics, Mar 2004; 20: 623 - 628.

p = inputParser;

addParameter(p, 'lambda', 20, @isnumeric)
addParameter(p, 'nbins', [], @isnumeric)
addParameter(p, 'plottype', 'scatter', @isstr)
addParameter(p, 'contourFlag', false, @islogical)
addParameter(p, 'msize', 5, @isnumeric)
addParameter(p, 'marker', 'o', @isstr)
addParameter(p, 'logx', false, @islogical)
addParameter(p, 'logy', false, @islogical)
addParameter(p, 'logc', false, @(x) islogical(x) || isnumeric(x))
addParameter(p, 'filled', true, @islogical)
addParameter(p, 'subsampling', 0, @isnumeric)
addParameter(p, 'maxF', NaN, @isnumeric)
addParameter(p, 'edges1', [], @isnumeric)
addParameter(p, 'edges2', [], @isnumeric)

parse(p,varargin{:});
p = p.Results;


if p.logx
    X = log10(X);
end
if p.logy
    Y = log10(Y);
end

minx = min(X,[],1);
maxx = max(X,[],1);
miny = min(Y,[],1);
maxy = max(Y,[],1);

if isempty(p.nbins)
    nbins = [min(numel(unique(X)),200) ,min(numel(unique(Y)),200) ];
else
    nbins = p.nbins;
end

if ~isempty(p.edges1)
    edges1 = p.edges1;
    ctrs1 = edges1(2:(end-1)) + .5*[diff(edges1(2:(end-1))) diff(edges1((end-2):(end-1)))];
else
    edges1 = linspace(minx, maxx, nbins(1)+1);
    ctrs1 = edges1(1:end-1) + .5*diff(edges1);
    edges1 = [-Inf edges1(2:end-1) Inf];
end

if ~isempty(p.edges2)
    edges2 = p.edges2;
    ctrs2 = edges2(2:(end-1)) + .5*[diff(edges2(2:(end-1))) diff(edges2((end-2):(end-1)))];
else
    edges2 = linspace(miny, maxy, nbins(2)+1);
    ctrs2 = edges2(1:end-1) + .5*diff(edges2);
    edges2 = [-Inf edges2(2:end-1) Inf];
end


[n,~] = size(X);
bin = zeros(n,2);
% Reverse the columns to put the first column of X along the horizontal
% axis, the second along the vertical.
[~,bin(:,2)] = histc(X,edges1);
[~,bin(:,1)] = histc(Y,edges2);
H = accumarray(bin,1,nbins([2 1])) ./ n;
G = smooth1D(H,nbins(2)/p.lambda);
F = smooth1D(G',nbins(1)/p.lambda)';
% = filter2D(H,lambda);

if p.logc
    F = F.^p.logc;
end
if isnan(p.maxF)
    maxF = max(F(:));
else
    maxF = p.maxF;
end
F = min(F/maxF,1);

if p.subsampling>0
    idx = randperm(length(X));
    idx = idx(1:min(p.subsampling,end));
else
    idx = 1:length(X);
end

switch lower(p.plottype)
    
    case 'surf'
        h = surf(ctrs1,ctrs2,F,'edgealpha',0);
    case 'mesh'
        h = mesh(ctrs1,ctrs2,F);
    case 'contour'
        [~, h] =contour(ctrs1,ctrs2,F);
    case 'image'
        nc = 256;
        colormap(repmat(linspace(1,0,nc)',1,3));
        h =image(ctrs1,ctrs2,floor(nc.*F) + 1);
    case 'scatter'
        ind = sub2ind(size(F),bin(:,1),bin(:,2));
        col = F(ind);
        if p.filled
            h = scatter(X(idx),Y(idx),p.msize,col(idx),p.marker,'filled');
        else
            h = scatter(X(idx),Y(idx),p.msize,col(idx),p.marker);
        end
    otherwise
        error('dscatter:UnknownPlotType',...
            'Unknown plot type: %s.',p.plottype);
        
end



if nargout > 0
    hAxes = get(h,'parent');
end
%%%% This method is quicker for symmetric data.
% function Z = filter2D(Y,bw)
% z = -1:(1/bw):1;
% k = .75 * (1 - z.^2);
% k = k ./ sum(k);
% Z = filter2(k'*k,Y);

function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;

