function DrugRates = Rates_from_CTM_AT(CTM, AT)
% rates: [ (static;toxic) X phase ] -> 0 for no effect, 1 for max effect

if exist('AT','var')
else
    AT = zeros(4);
end

DrugRates = [-diag(CTM)'; diag(AT)'];