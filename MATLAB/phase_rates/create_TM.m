function [CTM, AT] = create_TM(CellRates)
% CREATE_TM creates two transition matrix:  CTM that contains rates of
%transition in G1,S,G2,M (kG1,etc..) and AT that contains rates of apoptosis in
%G1,S,G2,M (taG1, etc..). The full transition matrix is TM = CTM - AT

%define matrix of state transitions between cell cycle phase
%  CTM(i,j) = from phase j to phase i:
    %   G1             S               G2              M
% CTM=[-kG1S          0               0               2*kMG1         ;     %G1
%      kG1S           -kSG2           0               0              ;     %S
% ....
CTM = -diag(CellRates(1,:)) + diag(CellRates(1,1:3),-1);
CTM(1,4) = 2*CellRates(1,4);

%Apoptotic matrix: defines the probability of apoptosis in a certain phase  
AT = diag(CellRates(2,:));

end

