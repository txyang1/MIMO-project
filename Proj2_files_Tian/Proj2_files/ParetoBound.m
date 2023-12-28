function [R] = ParetoBound(H,P,S)
% Function [R] = ParetoBound(H,P,S)
%
% The function calculates S coordinate points at the
% Pareto Boundary of the two-user MIMO MAC capacity region
% via varying the weights w_1 and w_2.
%
% Inputs
% H: M x N x 2 array of given channels to the users
% P: 2 x 1 column vector of available transmit power
% S: number of Pareto Boundary sample points S
% Outputs
% R: 2 x 2S matrix of rate region coordinates

M = size(H,1);
N = size(H,2);
