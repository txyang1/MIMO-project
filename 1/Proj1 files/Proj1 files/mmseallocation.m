function [psi,mu,K] = mmseallocation(phi,Ptx)
%
% MMSE minimizing power allocation procedure under the total
% transmit power constraint sum(psi)<=Ptx.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Ptx: available sum transmit power Ptx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% mu: value of the optimal waterlevel mu^bar
% K: number of active (non-zero) streams K

phi = phi(:);

phi_frac = 1./phi;
[phi_frac, index] = sort(phi_frac);
N = length(phi);

for K = 1:N
    mu = (Ptx + sum(phi_frac(1:K))) / sum(sqrt(phi_frac(1:K)));
    if K == N || mu < sqrt(phi_frac(K+1))
        break;
    end
end

psi = zeros(N,1);
for i = 1:K
    psi(i) = mu * sqrt(phi_frac(i)) - phi_frac(i);
end

pi = eye(N);
pi_permu = pi(index, :);
psi = pi_permu' * psi;

end
