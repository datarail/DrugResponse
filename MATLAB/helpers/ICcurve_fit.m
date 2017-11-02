%
% function [xI50, Hill, Einf, xMax, Area, r2, EC50, fit_final, p, log, flag] = ...
%     ICcurve_fit(Conc, Values, fit_type, opt)
%
%   Inputs:
%   ---------
%
%   Conc:   expected in uM, NOT logged
%   Growth: normalized to control (=1) fit_type
%   fit_type:   'GR50' (default if any Values<0), 'GI50' or 'IC50' (second default)
%
%   options:
%       - plotting: plot the curves [false]
%       - priors:   seed for the fitting
%       - ranges:   range for the fitting
%       - fitting:  average/individual for replicates [average]
%       - pcutoff:  cutoff for F-test [0.05]
%       - capped:   cap points with enhanced growth [true]
%       - extrapolrange: maximum fold to extrapolate xI50 [10]
%
%   Outputs:
%   ----------
%   concentrations are output in uM
%
%   IC50 => Growth is relative to control (end/ctrl)
%   GI50 => Growth is relative to growth of the control:
%               (end-day0)/(ctrl-day0)
%   GR50 => normalized growth is relative to growth of the control
%           Growth should not comprise the control value; each replicate is
%           a column
%
%   Area => sum of (1-cell count) devided by range --> average per order of
%   magnitude.
%
%   flag = 1 if fit is successful, 2 if Robust was used, 0 if linear used
%

function [xI50, Hill, Einf, xMax, Area, r2, EC50, fit_final, p, log, flag] = ...
    ICcurve_fit(Conc, Values, fit_type, opt)


% parameters : Einf  EC50 (uM)   HS 
priors = [.1 median(Conc) 2];

ranges = [
    0 1    %Einf
    max(min(Conc)*1e-4,1e-7) min(max(Conc)*1e2, 1e3)  %E50
    .1 5    % HS
    ]';
ranges(:,2) = -log10(ranges([2 1],2));
priors(2) = -log10(priors(2));

plotting = 0;
fitting = 'average';
pcutoff = .05;
capped = 1.05;
extrapolrange = 10;
Robust = [];

if ~exist('fit_type','var') || isempty(fit_type)
    if any(Values<0)
        fit_type = 'GR50';
    else
        fit_type = 'IC50';
    end
end

% 'IC50'
switch fit_type
    case 'IC50'
        ranges(1,1) = 0; % lowest Einf
    case 'GI50'
        ranges(1,1) = -2; % lowest GIinf; assuming that cells are at least 50% more than seeding.
    case 'GR50'
        ranges(1,1) = -1; % lowest GRinf;

end

% override default options
if exist('opt','var')
    fields = {'plotting', 'priors', 'ranges', 'fitting', 'pcutoff' 'capped' ...
        'extrapolrange' 'Robust'};
    for field = fields
        if isfield(opt,field{:})
            eval([field{:} ' = opt.' field{:} ';'])
        end
    end
end

Conc = ToRow(Conc);
if isvector(Values)
    Values=ToColumn(Values);
    fitting = 'average';
end
assert(all([1 size(Values,1)]==size(Conc)))


switch fitting
    case 'average'
        g = mean(Values,2)';
        g = g(sortidx(Conc));
        Conc = Conc(sortidx(Conc));
    case 'individual'
        for i=1:size(Values,2)
            [xI50(i), Hill(i), Einf(i), Area(i), r2(i), EC50(i), fit_final{i}, p(i), log{i}] = ...
                ICcurve_fit(Conc, Values(:,i)', fit_type, opt);
        end
        return

    otherwise
        error('wrong input for ''fitting''')

end



% remove the case of enhanced proliferation to avoid failure of F-test
if capped>0
    Conc = Conc(~isnan(g));
    g = g(~isnan(g));
    g = min(g, capped);
end

Npara = 3; % N of parameters in the growth curve
log = '';
flag = 1;
if isempty(Robust) || ~Robust
    [fit_res, gof] = sigmoidal_fit(Conc,g,'off');
else
    [fit_res, gof] = sigmoidal_fit(Conc,g,'off');
    r2 = gof.rsquare;
    if r2<.7 && length(g)>ceil(Npara) % don't use it is there are 4 points.
        warnprintf('Using robust fit (previous r=%.2f)', gof.rsquare)
        [fit_res, gof] = sigmoidal_fit(Conc,g,'Bisquare');
        log = sprintf('Robustfit (r=%.2f) -> ', gof.rsquare);
        fprintf('Robustfit (r=%.2f) -> ', gof.rsquare);
        flag = 2;
        if r2>gof.rsquare
            [fit_res, gof] = sigmoidal_fit(Conc,g,'off');
            fprintf('discared for normal fit\n');
        else
            fprintf('accepted \n');
        end
    end
end
[fit_res_flat, gof_flat] = flat_fit(Conc,g);

% F-test for the models
Npara_flat = 1;

RSS2 = gof.sse;
RSS1 = gof_flat.sse;

df1 = (Npara -Npara_flat);
df2 = (length(g) -Npara +1);
F = ( (RSS1-RSS2)/df1 )/( RSS2/df2 );
p = 1-fcdf(F, df1, df2);

xc = 10.^(log10(min(Conc)/extrapolrange):.05:log10(max(Conc)*extrapolrange));
r2 = gof.rsquare;

% results
Area = sum( (1-(g(2:end)+g(1:(end-1)))/2) .* diff(log10(Conc))) / ...
    diff(log10(Conc([1 end]))); % normalized version of the AUC
xMax = min(g(end-[1 0])); % robust minimum on the last 2 concentrations

if p>=pcutoff || isnan(RSS2) % failure of robust fit
    xI50 = +Inf;
    EC50 = +Inf;
    Hill = 0;
    Einf = min(g(end-[1 0])); % robust minimum on the last 2 concentrations
    log = [log '**** USING LINEAR FIT **** r2= ' num2str(gof_flat.rsquare,'%.2f')];
    fit_final = fit_res_flat;

    flag = 0;

else
    log = [log 'r2 = ' num2str(gof.rsquare,'%.2f')];

    fit_growth = fit_res(xc);

    Einf = fit_res.a;
    EC50 = 10^-fit_res.b;
    Hill = fit_res.c;
    fit_final = fit_res;

    xI50 = EC50*( ( ( (1-Einf)/(.5-Einf) )-1) ^(1/Hill));
    if any(fit_growth<.5) && any(fit_growth>.5) % inter/extrapolation is fine
        log = [log '\t' fit_type ' in the range of data'];
    elseif all(fit_growth>.5)
        if xI50>extrapolrange*max(Conc) || imag(xI50)~=0
            xI50 = Inf;
            log = [log '\t' fit_type '>' num2str(extrapolrange) '*max(Conc) --> +Inf'];
        end
    elseif all(fit_growth<.5)
        if xI50<min(Conc)/extrapolrange || imag(xI50)~=0
            xI50 = -Inf;
            log = [log '\t' fit_type '<min(Conc)/' num2str(extrapolrange) ' --> -Inf'];
        end
    else
        xI50 = NaN;
        warning(['undefined ' fit_type])
        log = [log '\tundefined ' fit_type ' --> NaN'];
    end

end


if plotting

    errorbar(log10(Conc), mean(Values,2), std(Values,[],2) ,'.k-');

    hold on
    plot(log10(xc), fit_res(xc),'-r')

    plot(log10(EC50)*[1 1], [0 .5+Einf/2], '.b-')
    plot(log10(Conc([1 end]))+[1 0], [1 1]*Einf, '.b-')

    plot(log10(xI50)*[1 1], [0 .5], '.b-')


    if p>=.05
        plot(log10(xc), fit_res_flat(xc),'-g');
    end

    title(sprintf('r^2 = %.3f', r2))
    xlim(log10([min(Conc)/extrapolrange extrapolrange*max(Conc)]))
    ylim([min(min(Values(:)), max(-.5,1-Einf)) max(Values(:))*1.05])
end


    function [fit_result, gof2] = sigmoidal_fit(doses, response, Robust)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',ranges(1,:),...
            'Upper',ranges(2,:),...
            'Startpoint',priors);
        % b' = -log10(b) -> search in the linear domain : (x/b) = (x*10^b')
        f = fittype('a + (1-a) ./ ( 1 + (x*(10.^b)).^c)','options',fitopt);
        [fit_result, gof2] = fit(doses', response',f,'Robust',Robust);
    end

    function [fit_result, gof2] = flat_fit(doses, response)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',ranges(1,1),...   % min Emax
            'Upper',1.2,...   % max E0
            'Startpoint',priors(1));
        f = fittype('a+0*x','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f);
    end

end
