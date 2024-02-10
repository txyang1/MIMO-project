function [psi,mu,K] = waterspilling(phi,E)
% Function [psi,mu,K] = waterspilling(phi,Etx)
%
% Computes the orthogonal projection of phi1,...,phiN onto
% a cone with constraints psik>=0 and sum(psi)<=E.
%
% Inputs
% phi: vector of eigenvalues phi1,...,phiN
% E: available sum transmit power E
% Outputs
% psi: vector of optimal power allocations psi1,...,psiN
% mu: value of the optimal waterlevel mu'
% K: number of non-zero data streams K

phi = phi(:);
phi(phi<0) = 0;

% determine N
N = length(phi);

% sort phi in descending order
phi_sort = sort(phi,'descend');

% compute the array mu_arr for all K=1:N
mu_arr = (cumsum(phi_sort)-E)./(1:N)';

% determine K
K = find(mu_arr < phi_sort,1,'last');

% compute mu and psi
mu = mu_arr(K);
psi = max(zeros(N,1), phi - mu);
