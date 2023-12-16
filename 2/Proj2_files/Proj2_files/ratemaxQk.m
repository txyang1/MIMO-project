function [Qk,Ck] = ratemaxQk(Xk,Pk)
%
% The function calculates the maximum rate and
% the rate maximizing transmit covariance.
%
% Inputs
% Xk: N x N matrix of effective noise covariance
% Pk: Maximum available transmit power
% Outputs
% Qk: single user rate maximizing transmit covariance matrix
% Ck: maximum single user rate

% TODO
N = size(Xk,1);
[Vk, Dk] = eig(Xk);
phi = real(diag(Dk));
[psi, ~, ~] = waterfilling(phi, Pk);
Qk = Vk * diag(psi) * Vk';
Ck = real(log2(det(eye(N) + Xk*Qk)));
