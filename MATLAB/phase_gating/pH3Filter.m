function [CCfrac, pH3CellIdentity, pH3cutoff, pH3lims, logpH3, logtxt] = pH3Filter(pH3, CellIdentity, varargin)
% [CCfrac, pH3CellIdentity, pH3cutoff, pH3lims] = pH3Filter(pH3, CellIdentity, ...)
%
%
% inputs are :
%  pH3  -> pH3 values
%  CellIdentity -> cell cycle identity (0=unclass, 1=G1, 2=S, 3=G2/M)
%  
% optional inputs are:
%  plotting     -> generates plots
%  interactive  -> prompt user to validate gating
%  savefigure   -> name to save image of the results
%
%  pH3cutoff    -> predefined cutoff
%
%  xpH3         -> sampling values for pH3 channel
%  pH3lims      -> plot range for pH3 channel
%
% outputs are:
%  CCfrac          -> fraction of cells in each phase (G1, S, G2, M, unclass)
%  pH3CellIdentity -> cell cycle identity (0=unclass, 1=G1, 2=S, 3=G2, 10.0/10.1/10.2/10.3=M)
%  pH3cutoff       -> selected cutoff for pH3 channel
%  pH3lims         -> selected range for pH3 channel
%
%

assert(all(size(pH3)==size(CellIdentity)))

p = inputParser;

addParameter(p, 'plotting', false, @islogical)
addParameter(p, 'interactive', false, @islogical)
addParameter(p, 'xpH3', 2.5:.02:8, @isvector)
addParameter(p, 'pH3lims', [], @(x) all(size(x)==[1 2]) && x(2)>x(1))
addParameter(p, 'pH3cutoff', [], @(x) isscalar(x) || isempty(x))
addParameter(p, 'savefigure', '', @ischar)

parse(p,varargin{:});
p = p.Results;

if p.interactive, p.plotting = true; end
if ~isempty(p.savefigure), p.plotting = true; end

%%


logpH3 = log10(min(max(pH3, 10^p.xpH3(3)),10^p.xpH3(end-2)));
% only take the cells in G1 or G2; to avoid contamination of pH3 by EdU
if ~any(CellIdentity==1 | CellIdentity==3)
    f = ksdensity(logpH3(CellIdentity==1 | CellIdentity==3), p.xpH3, 'width', ...
        4*diff(p.xpH3(1:2)));
else
    % % exceptional case ... take all cells
    f = ksdensity(logpH3, p.xpH3, 'width', ...
        4*diff(p.xpH3(1:2)));
end

if ~isempty(p.pH3cutoff)
    % already defined gate
    pH3cutoff = p.pH3cutoff;
    logtxt = 'pH3: used predefined cutoff';
end

if isempty(p.pH3cutoff) || mean(logpH3>=pH3cutoff)>.1 % if too many M phase cells switch to adaptive mode
    logtxt = 'pH3: adaptive pH3 cutoff';
    if ~isempty(p.pH3cutoff), logtxt = [logtxt ' (predefined has too many M-cells)']; end
    

    % determine the spread of the pH3 data to define cutoff
    [~, pk, pH3wdth] = findpeaks(f,'npeaks',3, ...
        'widthreference','halfprom','sortstr','descend');
    % enforce that no more than 30% of cells are in M-phase
    minidx = find((cumsum(f)/sum(f))>.3,1,'first')-5;
    if any(pk>=minidx)
        pH3wdth = pH3wdth(find(pk>=minidx,1,'first'));
        pk = max(pk(find(pk>=minidx,1,'first')), find((cumsum(f)/sum(f))>.3,1,'first'));
    else
        pk = minidx; pH3wdth = max(pH3wdth);
    end
    
    [~, minpk] = findpeaks(-f(pk:end),'npeaks',1, 'minpeakheight', -.6*max(f));
    minpk=minpk+pk-1;
    pH3cutoff = p.xpH3(ceil(max(min(minpk, pk+9*pH3wdth), pk+2*pH3wdth)));
    if isempty(pH3cutoff)
        pH3cutoff = p.xpH3(find(smooth(f,5)>.1/length(pH3),1,'last')+1); end
end

if ~isempty(p.pH3lims), pH3lims = p.pH3lims; else
    pH3lims = quantile(logpH3, [5e-3 .995])+[-1 10]*3*diff(p.xpH3(1:2)); end
if max(pH3lims)<pH3cutoff
    pH3lims(2) = pH3cutoff+.02; end

if p.plotting
    % plotting results
    
    plot_pos = [
    .07 .3 .4 .67;
    .6 .65 .3 .3
    .6 .08 .3 .5];

    get_newfigure(45679,[5 5 550 270])
    annotation('textbox', [.02 .02 .96 .08], 'string', logtxt)
    
    % plot original data
    get_newaxes(plot_pos(1,:),1)
    plot(p.xpH3, log10(f+max(f)/100)-log10(max(f)/100))
    
    fall = ksdensity(logpH3, p.xpH3, 'width', 2.5*diff(p.xpH3(1:2)));
    plot(p.xpH3, log10(fall+max(fall)/100)-log10(max(fall)/100), '--')
    
    plot([1 1]*pH3cutoff, ...
        [0 .5]*log10(max(f)), '-', 'color', [.6 .9 .1]);
    pltgt1 = plot([1 1]*pH3cutoff, ...
        [0 .5]*log10(max(f)), 'r-');
    if ~isempty(p.pH3cutoff) && p.pH3cutoff~=pH3cutoff
        plot([1 1]*p.pH3cutoff, ...
            [0 .5]*log10(max(f)), '-', 'color', [.6 .6 .6]);
    end
    
    xlim(pH3lims)
    ylim([0 log10(max(f))-log10(max(f)/100)+.1])
    set(gca,'xtick',[],'ytick',[])
    
    pieax = get_newaxes(plot_pos(2,:));   
    set(gca,'xtick',[],'ytick',[],'visible','off') 
    axis square
    
    pieax2 = get_newaxes(plot_pos(3,:));   
    set(gca,'xtick',[],'ytick',[],'visible','off') 
    axis square
end


%% finalize the results
EvalMphase();



if p.interactive
    figpos = get(gcf,'position');
    
    minpH3 = uicontrol('style', 'slider', 'callback', {@setCutoff,1});
    minpH3.Units = 'normalized';
    minpH3.InnerPosition = [plot_pos(1,1)-15/figpos(3) plot_pos(1,2)-.04 plot_pos(1,3)+30/figpos(3) .03];
    minpH3.Value = max(0,(pH3cutoff-pH3lims(1))/diff(pH3lims));
    
    
    
    approve = uicontrol('style', 'pushbutton');
    approve.Units = 'normalized';
    approve.Position = [.1 .03 .3 .13];
    approve.String = 'Approve';
    approve.Callback = @approveGate;
    
    waitfor(approve, 'backgroundcolor', 'g')
    
end

if p.plotting
    if ~isempty(p.savefigure)
        set(gcf,'Renderer','painters')
        saveas(gcf,p.savefigure)
    end
    
end
%%
    function setCutoff(src, event, x)
        pH3cutoff = (diff(pH3lims)*src.Value)+pH3lims(1);        
        set(pltgt1, 'XData', pH3cutoff*[1 1])
        logtxt = 'pH3: Manual adjustment';
        EvalMphase();
    end


    function Midx = EvalMphase()
        Midx = logpH3>=pH3cutoff;
        pH3CellIdentity = CellIdentity;   
        pH3CellIdentity(Midx) = 4+ pH3CellIdentity(Midx)/10;
        
        for id = 1:5
            CCfrac(id) = mean(floor(pH3CellIdentity)==mod(id,5));
        end
        
        if p.plotting
            set(gcf,'currentaxes', pieax)
            cla
            ptxt = pie(mean([Midx ~Midx])+1e-4 ,{sprintf('M-phase cells (%.1f%%)', ...
                100*mean(Midx)), 'Other'});
            set(ptxt(2),'fontsize',12, 'fontweight','bold')
            
            set(gcf,'currentaxes', pieax2)
            cla
            pie(CCfrac+1e-4, {'G1' 'S' 'G2' 'M' 'other'})
        end
    end

    function approveGate(src, event)
        minpH3.Visible = 'off';
        
        set(src, 'backgroundcolor', 'g')
    end


end