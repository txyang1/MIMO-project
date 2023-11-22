function [psi,mu,K] = powerMinimization(phi,r)
%
% Power Minimizing allocation under a minimum rate constraint
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% r: minimum rate requirement
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% mu: value of the optimal waterlevel mubar
% K: number of active (non-zero) data streams K
phi = phi(:);
phi= sort(phi,'descend');
phi_inv = 1./phi(:);
%find mu optimal and K
 for K=1:length(phi_inv)
     mu = ((2^r)/prod(phi(1:K)))^(1/K);
     if K==length(phi_inv)|| mu<phi_inv(K+1)
         break;
     end    
 end
%find psi
psi = zeros(size(phi));
for i = 1:K
    psi(i) = max(0, (mu - phi_inv(i)));
end


% TODO

