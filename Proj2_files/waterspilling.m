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
% determine N
N = length(phi);

% sort phi in descending order
phi_sort = sort(phi,'descend');

% compute the array mu_arr for all K=1:N
mu_arr = (cumsum(phi_sort)-Ptx)./(1:N)';

% determine K
K = find(phi_sort-mu_arr>0,1,'last');

% compute mu and psi
mu = mu_arr(K);
psi = max(zeros(N,1), phi_sort -mu_arr);

%Team members: Tingxin Yang, Tian Yu
