function [Qk,Ck] = ratemaxQk(Xk,Pk)
%
% The function calculates the maximum rate and
% the rate maximizing transmit covariance.
%
% Inputs
% Xk: N x N matrix of effective inverse noise covariance
% Pk: Maximum available transmit power
% Outputs
% Qk: single user rate maximizing transmit covariance matrix
% Ck: maximum single user rate

% TODO
N = size(Xk,1);

phi = real(eig(Xk)); %eigenvalues
[V,~] = eig(Xk); %modal matrix
[psi,~,~] = waterfilling(phi,Pk);
Qk = V*diag(psi)*V';
Ck = real(log2(det(eye(N)+Xk*Qk)));

