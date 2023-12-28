function [psi,mu,K] = waterspilling(phi,Ptx)
% function [psi,mu,K] = waterspilling(phi,Ptx)
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
phi(phi<0) = 0;

%TODO
phi= sort(phi,'descend');

for j = 1:numel(phi)
    s_temp(j) = max(phi(j),0);
end

if sum(s_temp)<=Ptx
    mu = 0;
else
    for K=1:length(phi)-1

        %Problem here with the index K?

        mu = 1/K*(sum(phi(1:K))-Ptx);
        if mu>=phi(K+1) && mu<phi(K)
            break;
        end
    end    
end

psi = zeros(size(phi));
for i = 1:numel(phi)
    psi(i) = max(0, (phi(i)-mu));
end

%Team members: Tingxin Yang, Tian Yu