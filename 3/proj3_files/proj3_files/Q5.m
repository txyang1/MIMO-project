clear all
close all
clc

load('exampleMIMOBC.mat');

order = [2,1];
S = MACtoBCtransform(Q,H,order);
[ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order );

% Ptx1q = trace(Q{1})
% Ptx2q = trace(Q{2})
% 
% Ptx1s = trace(S{1})
% Ptx2s = trace(S{2})
% 
% Rsum = sum(R_BC)