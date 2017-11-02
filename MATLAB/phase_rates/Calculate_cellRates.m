function CellRates = Calculate_cellRates(k0Rates, conc, MaxDrugRates, ECDrugRates)
% CREATE_TM creates two transition matrix:  CTM that contains rates of
%transition in G1,S,G2,M (kG1,etc..) and AT that contains rates of apoptosis in
%G1,S,G2,M (taG1, etc..). The full transition matrix is TM = CTM - AT

% fixed hill coefficient
n = 1.5;

% Hill function
DrugEffect = MaxDrugRates .* ((conc.^n)./((conc.^n)+(ECDrugRates.^n)));

CellRates = [k0Rates(1,:).*(1-DrugEffect(1,:))  % multiplicative for the static part
    DrugEffect(2,:)];           % additive for the toxic part
