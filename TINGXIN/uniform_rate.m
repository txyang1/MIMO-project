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
psi = zeros(size(phi));
N = length(phi);


  for K = 1:N
      if K < numel(phi) && sum(log(1 + phi(1:K) * Ptx / K)) > sum(log(1 + phi(1:K+1) * Ptx / (K+1)))
        break;
      end 
  end

 psi(1:K) = Ptx/K;
end



% TODO