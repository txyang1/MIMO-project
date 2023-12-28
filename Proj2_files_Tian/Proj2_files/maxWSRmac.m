function [Cwsum,Q] = maxWSRmac(H,P,w)
% Function [Cwsum,Q] = maxWSRmac(H,P,w)
%
% The function calculates the WSR capacity Cwsum and
% the optimal tranmit covariance matrices Q of the
% two-user MAC via a gradient projection method.
%
% Inputs
% H: M x N x K array of given channels to the users
% P: K x 1 column vector of available transmit power
% w: K x 1 column vector of weights w_1,...,w_K
% Outputs
% Cwsum: maximum WSR after convergence (WSR capacity)
% Q: N x N x K array of optimal covariance matrices

M = size(H,1);
N = size(H,2);
K = size(H,3);

fun = @(X) -wsrQ(X,H,w);
grad = @(X) wsrGradQ(X,H,w);
proj = @(X) projQ(X,P);

[Q,fval] = projGrad(fun,grad,proj,zeros(N,N,K),1e3);

Cwsum = -fval;
