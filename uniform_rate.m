function [psi,K] = uniform_rate(phi,Ptx)
%
% Rate maximizing uniform power allocation for the given average
% total transmit power constraint sum(psi)<=Ptx.
%
% Inputs
% phi: vector of eigenmode coefficients phi1,...,phiN
% Ptx: available sum transmit power Ptx
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% K: number of active data streams K

phi = phi(:);
phi = sort(phi,'descend');

N = numel(phi);

for K = 1:N
    if (K ~= N) && (sum(log2(1 + Ptx/K*phi(1:K))) <= sum(log2(1 + Ptx/(K+1)*phi(1:K+1))))
        continue
    else
        break
    end
end

psi = zeros(size(phi));
psi(1:K) = Ptx/K; %power allocation only till index K, from K+1 psi is 0

%Team members: Tian Yu, Tingxin Yang