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

% TODO
