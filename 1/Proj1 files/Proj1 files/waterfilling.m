function [psi,mu,K] = waterfilling(phi,Ptx)
% Function [psi,mu,K] = waterfilling(phi,Ptx)
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
phi_frac = 1./phi(:);
[phi_frac, index] = sort(phi_frac);

for K = 1:length(phi)
    mu = 1/K * (Ptx + sum(phi_frac(1:K)));
    if K == length(phi_frac) || mu < phi_frac(K+1)
        break;
    end
end

psi = zeros(size(phi));
for i = 1:K
    psi(i) = mu - phi_frac(i);
end

pi = eye(length(index));
pi_permu = pi(index, :);
psi = pi_permu' * psi;


end
% TODO
