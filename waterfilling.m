function [psi,mu,K] = waterfilling(phi,Ptx)
%
% Waterfilling power allocation procedure under the total transmit
% power constraint sum(psi)<=Ptx.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Ptx: available sum transmit power Ptx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% mu: value of the optimal waterlevel mu^prime
% K: number of active (non-zero) data streams K

phi = phi(:);
%phi= sort(phi,'descend');
phi_inv = 1./phi(:);
phi_inv = sort(phi_inv, 'descend');

%find mu optimal and K
 for K=1:length(phi_inv)
     mu = 1/K*(Ptx+sum(phi_inv(1:K)));
     if K==length(phi_inv)|| mu<phi_inv(K+1)
         %print(K);
         %print(mu);
         break;
     end    
 end
 
%find psi
psi = zeros(size(phi));

for i = 1:K
    psi(i) = max(0, (mu - phi_inv(i)));
end

%psi = diag(psi);
