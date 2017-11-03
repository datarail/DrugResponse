function [CCpeaks, CCfrac, DNAGates, EdUGates, CellIdentity, logDNA, logEdU, DNAlims, EdUlims, logtxt] = CCphases(DNA, EdU, varargin)
% [CCpeaks, CCfrac, DNAGates, EdUGates, CellIdentity, logDNA, logEdU, DNAlims, EdUlims] = CCphases(DNA, EdU, ...)
%
% inputs are :
%  DNA  -> Hoechst intensity values
%  EdU  -> EdU intensity values
%
% optional inputs are:
%  plotting     -> generates plots
%  interactive  -> prompt user to validate gating
%  savefigure   -> name to save image of the results
%
%  DNAGates     -> predefined DNA gates
%  EdUGates     -> predefined EdU gates
%  CCseeds      -> seeds for DNA/EdU peaks
%
%  xDNA         -> sampling values for DNA channel
%  xEdU         -> sampling values for EdU channel
%  DNAlims      -> plot range for DNA channel
%  EdUlims      -> plot range for EdU channel
%
% outputs are:
%  CCpeaks      -> location of DNA/EdU peaks (G1, S, G2) x (DNA, EdU)
%  CCfrac       -> fraction of cells in each phase (G1, S, G2, unclass)
%  DNAGates     -> selected DNA gates (3 values)
%  EdUGates     -> selected EdU gates (2 values)
%  CellIdentity -> cell cycle identity (0=unclass, 1=G1, 2=S 2.1=droppped S, 3=G2/M)
%  logDNA       -> transformed values for DNA channel
%  logEdU       -> transformed values for EdU channel
%  DNAlims      -> selected range for DNA channel
%  EdUlims      -> selected range for EdU channel
%
%

assert(all(size(DNA)==size(EdU)))

p = inputParser;

addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'interactive', false, @islogical)
addParameter(p, 'xDNA', 2.5:.02:8, @isvector)
addParameter(p, 'xEdU', -.2:.02:5.3, @isvector)
addParameter(p, 'nsmooth', 5, @isnumeric)
addParameter(p, 'CCseeds', [], @(x) ismatrix(x) && all(size(x)==[3 2]))
addParameter(p, 'DNAGates', [], @(x) ismatrix(x) & length(x)==4 & all(~isnan(x(:))))
addParameter(p, 'EdUGates', [], @(x) isvector(x) & length(x)==2 & all(~isnan(x(:))))
addParameter(p, 'savefigure', '', @ischar)
addParameter(p, 'EdUlims', [], @(x) all(size(x)==[1 2]) && x(2)>x(1))
addParameter(p, 'DNAlims', [], @(x) all(size(x)==[1 2]) && x(2)>x(1))
addParameter(p, 'DNAwidth', [], @(x) isscalar(x) || isempty(x))


parse(p,varargin{:});
p = p.Results;
if p.interactive, p.plotting = true; end
if ~isempty(p.savefigure), p.plotting = true; end


%% determine the spread of the EdU data to calibrate cutoffs
e = -200:4e3;
f = ksdensity(EdU, e);
% find the EdU peak for G1/G2
[~,pk,wdth]=findpeaks(f,'npeaks',2,'widthreference','halfprom',...
    'sortstr','descend','minpeakheight', max(f)/10);
wdth = wdth(argmin(pk));
pk = e(ceil(min(pk)));

if any(EdU>pk+30)
    f2 = ksdensity(EdU(EdU>pk+30), e);
    [~,m]=findpeaks(-f2,'npeaks',2,'widthreference','halfprom','sortstr','descend');
    m = e(ceil(m(argmin(abs(m-500)))));
    m = max([m; pk+3*wdth]);
else
    m = pk+3*wdth;
end

offsetEdU = max(pk-1.5*wdth,1);

% log10 domain and capping
logDNA = log10(min(max(DNA, 10^p.xDNA(3)),10^p.xDNA(end-2)));
logEdU = log10(min(max(EdU-offsetEdU, 10^p.xEdU(3)),10^p.xEdU(end-2)));

% expected maximum EdU value for G1 (in log10), at least 20% of cells
maxEdU = max(log10(m-offsetEdU), quantile(logEdU, .2)+.1);
% expected minumum EdU value for S (in log10)
minEdU = max(log10(pk+2*wdth-offsetEdU), maxEdU-.1);
% expected difference between G1/G2 and S (in log10)
EdUshift = max(log10(pk+2*wdth-offsetEdU)-log10(max(pk-offsetEdU,1)),1);


%
if p.plotting
    %     currfig = gcf;
    
    get_newfigure(45654,[5 285 550 600])
    % define positions
    plot_pos = [
        .04 .81  .2 .17
        .04 .6 .2 .16
        .29  .81  .3 .17
        .29  .6 .3 .16
        .63 .6 .35 .38
        .12  .13  .48  .4
        .7 .3  .25  .2
        .7  .15  .25 .08];
    ax = [];
    
    ax(7) = get_newaxes(plot_pos(7,:),1);
    axis square
    set(gca,'xtick',[],'ytick',[],'visible','off')
    
    % plot original data
    ax(2) = get_newaxes(plot_pos(2,:),1);
    idx = randperm(length(EdU),min(length(EdU),1000));
    plot(EdU(idx), logEdU(idx),'.c')
    hold on
    plot(e, max(p.xEdU)*f/max(f), 'k-')
    plot(e, max(p.xEdU)*f2/max(f2), 'k--')
    plot([-200 300 NaN -100 500 ], ...
        [minEdU*[1 1] NaN EdUshift*[1 1]+log10(max(pk-offsetEdU,1)) ], '-r')
    plot(offsetEdU*[1 1], [0 5], ':r')
    plot([-100 m m], [maxEdU*[1 1] 0], '--r')
    ylim(p.xEdU([1 end]))
    xlim([-200 max(pk+5*wdth, 500)])
    
    
    % plot original data
    ax(1) = get_newaxes(plot_pos(1,:),1);
    f = ksdensity(DNA, 0:100:1e6);
    plot(0:100:1e6, f, '-k')
    
    % plots with both channels
    ax(6) = get_newaxes(plot_pos(6,:),1);
    dscatter(logDNA, logEdU, 'MSIZE', 20, 'marker', 'o')
    hold on
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first pass at the finding peaks by combining DNA and EdU channels



xDNA2 = p.xDNA;
xEdU2 = p.xEdU;
nbins = [length(xDNA2) length(xEdU2)]-1;
bin = NaN(length(logDNA),2);
[~,bin(:,2)] = histc(logDNA,xDNA2);
[~,bin(:,1)] = histc(logEdU,xEdU2);

% counting and smoothing
H = accumarray(bin,1,nbins([2 1]));
H(H==1) = 0;
H = H/length(logDNA);
G = smooth1D(H,p.nsmooth);
F = smooth1D(G',p.nsmooth)';

% finding peaks
Pk2D = imregionalmax(F);
[x,y] = find(Pk2D);

% store candidate peaks
prePksCandidates = [xDNA2(y)' xEdU2(x)' F(Pk2D)];
PksCandidates = flipud(sortrows(...
    prePksCandidates( (prePksCandidates(:,2)>(p.nsmooth+2)*diff(xEdU2([1 2]))) & ...
    prePksCandidates(:,3)>1e-5 & prePksCandidates(:,3)>max(prePksCandidates(:,3)/30),:),3));
% filter out the small peaks and one with EdU=0; sort descending

% less smoothing if only one peak found
if size(PksCandidates,1)<2
    G = smooth1D(H,p.nsmooth/2);
    F = smooth1D(G',p.nsmooth/2)';
    
    % finding peaks
    Pk2D = imregionalmax(F);
    [x,y] = find(Pk2D);
    
    % store candidate peaks
    prePksCandidates = [xDNA2(y)' xEdU2(x)' F(Pk2D)];
    PksCandidates = flipud(sortrows(...
        prePksCandidates( (prePksCandidates(:,2)>(p.nsmooth+2)*diff(xEdU2([1 2]))) & ...
        prePksCandidates(:,3)>1e-5 & prePksCandidates(:,3)>max(prePksCandidates(:,3)/20),:), 3));
    % filter out the small peaks and one with EdU=0; sort descending
    
elseif size(prePksCandidates,1)>2*size(PksCandidates,1) || size(prePksCandidates,1)>5
    % if too many peaks discarded -> more smoothing
    G = smooth1D(H,p.nsmooth*2);
    F = smooth1D(G',p.nsmooth*2)';
    
    % finding peaks
    Pk2D = imregionalmax(F);
    [x,y] = find(Pk2D);
    
    % store candidate peaks
    prePksCandidates = [xDNA2(y)' xEdU2(x)' F(Pk2D)];
    PksCandidates = flipud(sortrows(...
        prePksCandidates( (prePksCandidates(:,2)>(p.nsmooth+2)*diff(xEdU2([1 2]))) & ...
        prePksCandidates(:,3)>1e-5 & prePksCandidates(:,3)>max(prePksCandidates(:,3)/60),:), 3));
    % filter out the small peaks and one with EdU=0; sort descending
end

if p.plotting
    % plotting the results
    ax(5) = get_newaxes(plot_pos(5,:),1);
    imagesc(xDNA2, xEdU2, F)
    scatter(PksCandidates(:,1), PksCandidates(:,2), 20+sqrt(PksCandidates(:,3))*100, 'ok')
    set(gca,'ydir','normal')
    xlim(quantile(xDNA2,[.25 .75]))
    ylim(quantile(xEdU2,[.1 .75]))
end

PhasesCandidates = NaN(3,2);

% assign the candidate peaks to G1, S, or G2
if any(PksCandidates(:,2)-min(PksCandidates(:,2))>EdUshift & ...
        PksCandidates(:,2)>(minEdU+.2*EdUshift))
    % there is enough differences on EdU to expect a S phase
    
    % get the S phase peak as the one with high EdU and high density within 1.5*log10(2)
    PhasesCandidates(2,:) = PksCandidates(argmax(...
        (PksCandidates(:,2)-min(PksCandidates(:,2))>EdUshift & PksCandidates(:,2)>(minEdU+.2*EdUshift))...
        .*(sum((abs(repmat(PksCandidates(:,1), 1, size(PksCandidates,1))-...
        repmat(PksCandidates(:,1), 1, size(PksCandidates,1))')<(log10(2)*.75)).*...
        (repmat(PksCandidates(:,3),1,size(PksCandidates,1))),1)'+PksCandidates(:,3))),[1 2]);
    
    % find the G1 peak now
    temp = PksCandidates(PksCandidates(:,1)<(PhasesCandidates(2,1)+log10(2)*.05) & ...
        PksCandidates(:,1)>(PhasesCandidates(2,1)-.75*log10(2)) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, :);
    if ~isempty(temp) % there is a likely G1 peak
        % assign the G1 peak as the one close to the S peak on the DNA axis
        % with high density
        PhasesCandidates(1,:) = temp(argmax(temp(:,3)),[1 2]);
    end
    % assign the G2 peak as the one closest to the S peak on the DNA axis
    temp = PksCandidates(PksCandidates(:,1)>nanmean(PhasesCandidates(1:2,1)) & ...
        PksCandidates(:,1)<nanmean(PhasesCandidates(1:2,1))+log10(2) & ...
        PksCandidates(:,2)<PhasesCandidates(2,2)-EdUshift, :);
    if ~isempty(temp) % there is a likely G2 peak
        PhasesCandidates(3,:) = temp(argmax(temp(:,3)),[1 2]);
    end
    
else
    % most likely no S peak
    % -> take the two highest peak and assign as G1 and G2 based on DNA
    if ~isempty(PksCandidates)
        if size(PksCandidates,1)==1
            % if only one peak, just assign it to G1
            PhasesCandidates(1,:) = PksCandidates(1,[1 2]);
        else
            % if 2 peaks, take the top ones that are enough apart
            Pksdist = (dist([PksCandidates(:,1) zeros(size(PksCandidates,1),1)]')>.6*log10(2)).* ...
                (repmat(PksCandidates(:,3),1,size(PksCandidates,1))+ ...
                repmat(PksCandidates(:,3)',size(PksCandidates,1),1));
            [PksIdx1, PksIdx2] = find(Pksdist==max(Pksdist(:)),1,'first');
            if PksCandidates(PksIdx1,1)>PksCandidates(PksIdx2,1)
                PhasesCandidates(1,:) = PksCandidates(PksIdx2,[1 2]);
                PhasesCandidates(3,:) = PksCandidates(PksIdx1,[1 2]);
            else
                PhasesCandidates(1,:) = PksCandidates(PksIdx1,[1 2]);
                PhasesCandidates(3,:) = PksCandidates(PksIdx2,[1 2]);
            end
        end
    end
end

if p.plotting
    % plotting the results
    phases = {'G1' 'S' 'G2'};
    for i=1:3
        text(PhasesCandidates(i,1), PhasesCandidates(i,2), phases{i}, ...
            'fontsize', 14, 'fontweight', 'bold', 'color', [.7 .2 .2], ...
            'horizontalalign','center')
    end
end

EvaluatedPhasesCandidates = PhasesCandidates;

if ~isempty(p.CCseeds)
    % if seeds for the peaks were provided, assign them
    PhasesCandidates(all(~isnan(p.CCseeds),2),:) = p.CCseeds(all(~isnan(p.CCseeds),2),:);
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now working with each channel sequentially


% first working with the DNA content
f = ksdensity(logDNA,p.xDNA);

if p.plotting
    ax(3) = get_newaxes(plot_pos(3,:),1);
    plot(p.xDNA, f)
end

% take only the cells with low EdU
f = ksdensity(logDNA(logEdU<minEdU+.2*EdUshift & logEdU<maxEdU),p.xDNA);
if p.plotting, plot(p.xDNA, f, '--'), end

[pks, idx] = findpeaks(f, 'sortstr', 'descend');
idx = idx(pks>max(pks/10)); % remove lesser peaks
DNAPks = p.xDNA(idx(1:min(3,length(idx)))); % take the 3 highest peaks with low EdU

% find the DNA peak for G1
if length(DNAPks)>1
    % more than one candidate
    if ~isnan(PhasesCandidates(1,1))
        % 2D analysis matching G1 peak
        DNAPks = DNAPks(argmin(abs(DNAPks-PhasesCandidates(1,1))));
    elseif ~isnan(PhasesCandidates(2,1))
        % 2D analysis matching S peak
        DNAPks = max(DNAPks(DNAPks<PhasesCandidates(2,1)));
    else
        % take the peak with lowest DNA value (most likely case in doubt)
        DNAPks = min(DNAPks);
    end
end

if isempty(DNAPks)
    % missing G1 peak -> get something by default
    DNAPks = nanmin(PhasesCandidates(:,1))-log10(1.2);
end
if p.plotting, plot(DNAPks, .1, 'xk'); end

%% %%%%%%%%%%%%%%%%%%%%
% now working with EdU
f = ksdensity(logEdU,p.xEdU);
if p.plotting
    ax(4) = get_newaxes(plot_pos(4,:),1);
    plot(p.xEdU, f)
end
% find the low EdU (G1 and early S)
lE = ( (logDNA>DNAPks-1) & (logDNA<DNAPks+.1) ) & ...
    logEdU>2*p.nsmooth*diff(p.xEdU(1:2)) & logEdU<maxEdU;
if ~any(lE)
    % relax constrain on cell with EdU below thereshold
lE = ( (logDNA>DNAPks-1) & (logDNA<DNAPks+.1) ) & logEdU<maxEdU;
end
f = ksdensity(logEdU(lE),p.xEdU);
N = histcounts(logEdU(lE),p.xEdU);
f([true; smooth(N,3)<=1/3]) = 0; % remove single cells
if p.plotting, plot(p.xEdU, f, '--'), end

[EdUint, idx] = findpeaks(smooth(f,p.nsmooth),'sortstr','descend','npeaks',2);

if ~isempty(idx)
    % take the highest peak or lowest EdU level (G1 + G2)
    idx = idx(EdUint>.3*max(EdUint));
    EdUPks = p.xEdU(idx(argmin(idx)));
else
    EdUPks = median(logEdU(lE));
end
%
% get the peak with high EdU values (S)
hE = (logDNA>DNAPks-log10(2)/2) & (logDNA<DNAPks+log10(2)*1.5) & logEdU>EdUPks+EdUshift*.8;
if any(hE)
    % found some cells likely in S phase
    f = ksdensity(logEdU(hE),p.xEdU);
    if p.plotting, plot(p.xEdU, f, ':'), end
    
    [pks, idx] = findpeaks(smooth(f,p.nsmooth),'sortstr','descend');
    hEdUPks = idx(pks>max(pks/10)); % remove lesser peaks
    
    if any(p.xEdU(hEdUPks)>(EdUPks+EdUshift))
        % should be at least EdUshift above the EdU value in G1 peak
        EdUPks = [EdUPks ...
            p.xEdU(hEdUPks(find(p.xEdU(hEdUPks)>(EdUPks+EdUshift),1,'first'))) EdUPks];
    else
        % set by default if no good candidate
        EdUPks = [EdUPks (EdUPks+EdUshift) EdUPks];
    end
    
    % use the distribution to find the minimum of EdU as cutoff (between peaks)
    f = ksdensity(logEdU,p.xEdU);
    [~,idx] = findpeaks(-f,'sortstr','descend');
    EdUcutoff = p.xEdU(idx(find(p.xEdU(idx)>EdUPks(1) & p.xEdU(idx)<EdUPks(2),1,'first')));
    if isempty(EdUcutoff)
        EdUcutoff = p.xEdU(argmin(smooth(f',p.nsmooth)' + ...
            (p.xEdU<EdUPks(1) | p.xEdU>EdUPks(2))));
    end
    EdUlims = [p.xEdU(3) min(EdUPks(2)+(EdUPks(2)-EdUcutoff), p.xEdU(end-1))];
    
    if p.plotting
        % plot the location of the peaks
        plot(EdUPks, .1, 'xk');
        plot(EdUcutoff, .1, 'xk');
    end
    
    % get back to the DNA to find the S location
    hE = (logDNA>DNAPks-log10(2)/2) & (logDNA<DNAPks+log10(2)*1.5) & logEdU>EdUcutoff;
    % get the distribution of cells in S phase
    if any(hE)
        f = ksdensity(logDNA(hE),p.xDNA);
        
        if p.plotting
            set(gcf,'currentaxes',ax(3))
            plot(p.xDNA, f, '-.')
        end
        
        [~, idx] = findpeaks(smooth(f,3*p.nsmooth),'sortstr','descend', 'Npeaks', 1);
        DNAPks = [DNAPks p.xDNA(idx)];
    else
        % set a default value
        DNAPks = [DNAPks DNAPks+.5*log10(2)];
    end
    
else
    % no cells found in S phase
    
    % set arbitrary value for high S
    EdUPks = [EdUPks (EdUPks+EdUshift) EdUPks];
    EdUcutoff = mean(EdUPks([1 2]));
    EdUlims = [-.02 min(EdUPks(2)+(EdUPks(2)-EdUcutoff)+.1, p.xEdU(end-1))];
    
    % set default DNA content (1.4 fold)
    DNAPks = DNAPks+[0 log10(2)*.5];
end

% at least width of 1 for the S phase box
EdUGates = [EdUcutoff EdUPks(2)+max(EdUPks(2)-EdUcutoff,1)];

if ~isempty(p.EdUlims), EdUlims = p.EdUlims;
else, EdUlims(2) = max(EdUlims(2), EdUGates(2)+.01); end
    

%% %%%%%%%%%%%%%%%%%%%%
% working with DNA again to find the G2

hD = logDNA>DNAPks(1)+.4*log10(2) & logEdU<EdUGates(1);
% get the cells likely to be in G2 phase

if any(hD)
    % found some cells likely in G2 phase
    f = ksdensity(logDNA(hD),p.xDNA);
    
    if p.plotting
        plot(p.xDNA, f, ':')
    end
    
    [pks, idx] = findpeaks(smooth(f,p.nsmooth), 'sortstr', 'descend');
    idx = idx(pks>max(pks/10)); % remove lesser peaks
    hDNAPks = p.xDNA(idx);
    % should be around log10(2) above the DNA peak in G1
    hDNAPks = hDNAPks(hDNAPks>(DNAPks(1)+.5*log10(2)));
    
    if length(hDNAPks)>1
        % more than one candidate
        if ~isnan(PhasesCandidates(3,1))
            % 2D analysis matching G2 peak
            hDNAPks = hDNAPks(argmin(abs(hDNAPks-PhasesCandidates(3,1))));
        else
            if ~isnan(PhasesCandidates(2,1)) && any(hDNAPks>PhasesCandidates(2,1))
                % 2D analysis matching S peak
                hDNAPks = hDNAPks(hDNAPks>PhasesCandidates(2,1));
            end
            % take the peak closest to a 2-fold(most likely case in doubt)
            hDNAPks = hDNAPks(argmin(abs(hDNAPks-DNAPks(1)-log10(2))));
        end
    end
    if isempty(hDNAPks)
        % no good candidate for G2; set by default
        DNAPks = [DNAPks DNAPks(1)+log10(2)];
    else
        DNAPks = [DNAPks hDNAPks];
    end
    
    % find the split between G1 and G2
    f = ksdensity(logDNA(logEdU<EdUGates(1)),p.xDNA);
    [~, idx] = findpeaks(-smooth(f,p.nsmooth), 'sortstr', 'descend');
    DNAcutoff = p.xDNA(idx(p.xDNA(idx)>DNAPks(1) & p.xDNA(idx)<DNAPks(3)));
    
    if isempty(DNAcutoff)
        % no good candidate; set by default
        DNAcutoff = min(max(DNAPks(2),DNAPks(1)+.02), DNAPks(3)-.02);
    elseif length(DNAcutoff)>1
        DNAcutoff = DNAcutoff(1);
    end
    
else
    % no cells found in G2 phase
    
    % set default values (2-fold for G2)
    DNAcutoff = DNAPks(1)+.3*log10(2);
    DNAPks = [DNAPks DNAPks(1)+log10(2)];
    
end

%% find the cells dropping in S-phase

% G1 width (take only the centered cells to estimate width)
hG1 = abs(logDNA-DNAPks(1))<.3*log10(2) & logEdU<EdUGates(1);
if sum(hG1)>10
    NormFitG1 = fitdist(logDNA(hG1), 'Normal');
    G1minwidth = max(NormFitG1.icdf(.9)-NormFitG1.icdf(.1),.05);
    G1lim = min(NormFitG1.icdf(.99), DNAcutoff-.1*log10(2));
else
    G1lim = min(DNAcutoff-.1*log10(2), mean([DNAcutoff, DNAPks(1)]));
    G1minwidth = .05;
end

hG2 = abs(logDNA-DNAPks(3))<.3*log10(2) & logEdU<EdUGates(1);
if sum(hG2)>10
    NormFitG2 = fitdist(logDNA(hG2), 'Normal');
    G2minwidth = max(NormFitG2.icdf(.9)-NormFitG2.icdf(.1),.05);
    G2lim = max(NormFitG2.icdf(.01), DNAcutoff+.1*log10(2));
else
    % case when too few cells
    G2lim = max(DNAcutoff+.1*log10(2), mean([DNAcutoff, DNAPks(3)]));
    G2minwidth = .05;
end
%% define the DNA range and gates
d1 = DNAcutoff-DNAPks(1);
d2 = DNAPks(3)-DNAcutoff;
if ~isempty(p.DNAlims), DNAlims = p.DNAlims; else
    DNAlims = [max(DNAPks(1)-3*d1, p.xDNA(2)) min(DNAPks(3)+3*d2, p.xDNA(end-1))]; end
if p.plotting
    plot(DNAPks, 0, 'xk');
    plot(DNAcutoff, 0, 'xk');
end


% defines the gates
DNAGates = [DNAPks(1)-d1 G1lim G2lim DNAPks(3)+d2];
% inflate if too small
if diff(DNAGates([1 2]))<G1minwidth
    DNAGates([1 2]) = mean(DNAGates([1 2]))+[-.6 .4]*G1minwidth;
end
if diff(DNAGates([3 4]))<G1minwidth
    DNAGates([3 4]) = mean(DNAGates([3 4]))+[-.4 .6]*G2minwidth;
end
if DNAGates(3)<DNAGates(2)
    DNAGates(2:3) = mean(DNAGates(2:3));
end

DNAlims = [min([DNAlims DNAGates-.1]) max([DNAlims DNAGates+.1])];

% compare with input in present
EvaluatedDNAGates = DNAGates;
if ~isempty(p.DNAGates)
    if nanmax(DNAGates - p.DNAGates)>.2
        warnprintf('Large differences between input DNA gates and calculates ones')
    end
    DNAGates = p.DNAGates;
end

EvaluatedEdUGates = EdUGates;
if ~isempty(p.EdUGates)
    if nanmax(EdUGates - p.EdUGates)>.2
        warnprintf('Large differences between input DNA gates and calculates ones')
    end
    EdUGates = p.EdUGates;
end

%% get the fraction of cells in each phase
[CCfrac, CCpeaks, CellIdentity] = EvalCC();
logtxt = 'DNA/EdU: automatic gating';

%%
if p.plotting
    % plot the evaluated phase positions
    phases = {'G1' 'S' 'G2'};
    for i=5:6
        for j=1:3
            set(gcf,'currentaxes', ax(i))
            text(EvaluatedPhasesCandidates(j,1), EvaluatedPhasesCandidates(j,2), phases{j}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', [.6 .9 .1], ...
                'horizontalalign','center')
            text(PhasesCandidates(j,1), PhasesCandidates(j,2), phases{j}, ...
                'fontsize', 14, 'fontweight', 'bold', 'color', [.1 .1 .1], ...
                'horizontalalign','center')
        end
    end
    % plot the evaluated gates
    plot(ax(1), 10.^EvaluatedDNAGates([1 1 1 2 2 2 3 3 3 4 4]), ...
        max(ylim(ax(1)))*[0 1 NaN 0 1 NaN 0 1 NaN 0 1], '--', 'color', [.6 .9 .1])
    plot(ax(3), EvaluatedDNAGates([1 1 1 2 2 2 3 3 3 4 4]), ...
        max(ylim(ax(3)))*[0 1 NaN 0 1 NaN 0 1 NaN 0 1], '--', 'color', [.6 .9 .1])
    plot(ax(4), EvaluatedEdUGates([1 1 1 2 2]), ...
        max(ylim(ax(4)))*[0 1 NaN 0 1], '--', 'color', [.6 .9 .1])
    for i=5:6
        plot(ax(i), EvaluatedDNAGates([1 1 4 4  1  1 2 2  1  3 3 4]), ...
            [-1 EvaluatedEdUGates([2 2]) -1 NaN EvaluatedEdUGates([1 1]) -1 ...
            NaN -1 EvaluatedEdUGates([1 1])], '--', 'color', [.6 .9 .1], 'linewidth', 2)
    end
    
    hgates = [];
    hphases = [];
    refreshGates()
    
    plot(ax(5),DNAPks, EdUPks, 'xk')
    
    xlim(ax(3),DNAlims)
    xlim(ax(4),EdUlims)
    for i=[5 6]
        xlim(ax(i),DNAlims)
        ylim(ax(i),EdUlims)
    end
    
    if p.interactive
        % Manual adjustments for DNA content
        figpos = get(gcf,'position');
        DNAslides = {};
        for i=1:4
            DNAslides{i} = uicontrol('style', 'slider', 'callback', {@setDNAGates,i});
            DNAslides{i}.Units = 'normalized';
            DNAslides{i}.Position = [plot_pos(6,1)-15/figpos(3) plot_pos(6,2)-.03*i-.005 plot_pos(6,3)+30/figpos(3) .028];
            DNAslides{i}.Value = (DNAGates(i)-DNAlims(1))/diff(DNAlims);
        end
        
        % Manual adjustments for EdU content
        figpos = get(gcf,'position');
        
        EdUslides = {};
        for i=1:2
            EdUslides{i} = uicontrol('style', 'slider', 'callback', {@setEdU, i});
            EdUslides{i}.Units = 'normalized';
            EdUslides{i}.Position = [plot_pos(6,1)-.04*i-.02 plot_pos(6,2)-15/figpos(4) .04 plot_pos(6,4)+30/figpos(4)];
            EdUslides{i}.Value = (EdUGates(i)-EdUlims(1))/diff(EdUlims);
        end
        
        approve = uicontrol('style', 'pushbutton');
        approve.Units = 'normalized';
        approve.Position = plot_pos(8,:);
        approve.String = 'Approve';
        approve.Callback = @approveGate;
        
        waitfor(approve, 'backgroundcolor', 'g')
    end
    
    if ~isempty(p.savefigure)
        set(gcf,'Renderer','painters')
        saveas(gcf,p.savefigure)
    end
    
    %     figure(currfig)
end


%%

    function setDNAGates(src, event, x)
        logtxt = 'DNA/EdU: Manual adjustment';
        DNAGates(x) = (diff(DNAlims)*src.Value)+DNAlims(1);
        % check for proper ordering
        if x<4 && DNAGates(x)>DNAGates(x+1)
            DNAGates(x) = DNAGates(x+1);
        elseif x>1 && DNAGates(x)<DNAGates(x-1)
            DNAGates(x) = DNAGates(x-1);
        end
        for iG=1:4
            DNAslides{iG}.Value = (DNAGates(iG)-DNAlims(1))/diff(DNAlims);
        end
        % re-evluate assignments
        [CCfrac, CCpeaks, CellIdentity] = EvalCC();
        if p.plotting, refreshGates(), end
    end


    function setEdU(src, event, x)
        logtxt = 'DNA/EdU: Manual adjustment';
        EdUGates(x) = (diff(EdUlims)*src.Value)+EdUlims(1);
        % check for proper ordering
        if EdUGates(1)>EdUGates(2)
            EdUGates(1)=EdUGates(2);
            EdUslides{1}.Value = (EdUGates(1)-EdUlims(1))/diff(EdUlims);
        end
        % re-evluate assignments
        [CCfrac, CCpeaks, CellIdentity] = EvalCC();
        if p.plotting, refreshGates(), end
    end


    function refreshGates()
        if isempty(hphases)
            for iA = 1:2
                set(gcf,'currentaxes', ax(iA+4))
                for iG=1:3
                    hphases(iA,iG) = text(DNAPks(iG), EdUPks(iG), phases{iG}, ...
                        'fontsize', 16, 'fontweight', 'bold', 'color', [.7 .2 0], ...
                        'horizontalalign','center');
                end
            end
        end
        if isempty(hgates)
            for iG=[1 3:6]
                hgates(iG) = plot(ax(iG),NaN, NaN, '-r','linewidth', 1+(iG>4));
            end
        end
        
        for iG=1:3
            for iA = 1:2
                set(hphases(iA,iG), 'position', [DNAPks(iG), EdUPks(iG)])
            end
        end
        
        set(hgates(1), 'xdata', 10.^DNAGates([1 1 1 2 2 2 3 3 3 4 4]), ...
            'ydata', max(ylim(ax(1)))*[0 1 NaN 0 1 NaN 0 1 NaN 0 1])
        set(hgates(3), 'xdata', DNAGates([1 1 1 2 2 2 3 3 3 4 4]), ...
            'ydata', max(ylim(ax(3)))*[0 1 NaN 0 1 NaN 0 1 NaN 0 1])
        set(hgates(4), 'xdata', EdUGates([1 1 1 2 2]), ...
            'ydata', max(ylim(ax(4)))*[0 1 NaN 0 1])
        for iG=5:6
            set(hgates(iG), 'xdata', DNAGates([1 1 4 4  1  1 2 2  1  3 3 4]), ...
            'ydata', [-1 EdUGates([2 2]) -1 NaN EdUGates([1 1]) -1 ...
            NaN -1 EdUGates([1 1])])
    
        end
        set(gcf,'currentaxes',ax(7))
        cla
        ptxt = pie(CCfrac+1e-4, {'G1' 'S' 'G2' sprintf('other %.0f%%', 100*CCfrac(4))});
        set(ptxt(end),'fontsize',12, 'fontweight','bold')
        
    end

    function [frac, pks, cellID] = EvalCC()
        
        cellID = (logDNA>=DNAGates(1) & logDNA<DNAGates(2) & logEdU<EdUGates(1)) ... % G1
            +2*(logDNA>=DNAGates(1) & logDNA<DNAGates(4) & logEdU>=EdUGates(1) & logEdU<EdUGates(2)) ... % S
            +2.1*(logDNA>=DNAGates(2) & logDNA<DNAGates(3) & logEdU<EdUGates(1)) ... % dropped S
            +3*(logDNA>=DNAGates(3) & logDNA<DNAGates(4) & logEdU<EdUGates(1)); % G2
        
        for id = 1:4
            frac(id) = mean(floor(cellID)==mod(id,4));
        end
        
        for iG=1:3
            if sum(cellID==iG)>10
                [~, pkidx] = findpeaks(smooth( ksdensity(logDNA(cellID==iG),p.xDNA), ...
                    3*p.nsmooth), 'sortstr','descend', 'Npeaks', 1);
                DNAPks(iG) = p.xDNA(pkidx);
                [~, pkidx] = findpeaks(smooth(ksdensity( logEdU(cellID==iG),p.xEdU), ...
                    3*p.nsmooth), 'sortstr','descend', 'Npeaks', 1);
                EdUPks(iG) = p.xEdU(pkidx);
            else
                DNAPks(iG) = mean(DNAGates([iG 2]));
                EdUPks(iG) = mean([EdUGates(1) (iG==2)*EdUGates(2)]);
            end
            EdUPks(2) = max(EdUPks(2), EdUGates(1)+.1);
        end
        
        pks = [DNAPks' EdUPks'];
    end

    function approveGate(src, event)
        for iB=1:4
            DNAslides{iB}.Visible = 'off';
        end
        for iB=1:2
            EdUslides{iB}.Visible = 'off';
        end
        
        set(src, 'backgroundcolor', 'g')
    end

end