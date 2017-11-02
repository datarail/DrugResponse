function [ss_growth_rate, ss_dist_cell_cycle] = calculate_steady_state(CTM, AT)
% CALCULATE_STEADY_STATE calculates the steady state distribution of a
% number of cells (ncells cycling according to rates in CTM (state transition matrix) and
% AT (apoptotic matrix)

if exist('AT','var')
    TM = CTM - AT;
else
    TM = CTM;
end
    

% calculate the steady state growth rate and cell fraction distribution
% using eigenvalues (TM * x = mu * x). The positive eigenvalue is the growth
% rate, the corresponding eigenvector is the fractional distribution in the
% cell cycle
[V,D] = eig(TM); %V: eigenvectors (columns0, D: eigenvalues (diagonal)
E=diag(D);
%find real eigenvalues
E_real=(imag(E)==0); 
%find positive eigenvalues
E_pos=(real(E)>=0); 
%find real eigenvectors
V_real=all(imag(V)==0)'; 
%find eigenvectors with coherent signs
V_sign=sign(real(V));
V_sign_pos=V_sign; V_sign_neg=V_sign;
V_sign_pos(V_sign_pos==0)=1;  V_sign_pos(V_sign_pos==-1)=0; 
V_sign_neg(V_sign_neg==0)=-1; V_sign_neg(V_sign_neg==1)=0; 
V_coherent=(all(V_sign_pos)|all(V_sign_neg))';
% if there is a unique real eigenvalue with a coherent eigenvector (or all
% are the same value), then that is the stable growth rate
i_sol=intersect(find(E_real==1),find(V_coherent));
if ((~isempty(i_sol) && (length(i_sol)==1)) || (sum(E==E(1))==length(E)))
    ss_growth_rate=E(i_sol(1));
    ss_dist_cell_cycle=abs(V(:,i_sol(1)))./sum(abs(V(:,i_sol(1)))); 
else
    % there isn't a single real eigenvalue with coherent eigenvector (either
    % none or multiple with different values), no stable growth rate
    ss_growth_rate=NaN;
    ss_dist_cell_cycle=repmat(NaN,size(TM,1),1);
end   



