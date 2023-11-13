clear all;
clc
load('example_channels.mat');
S = H' * inv(Cn) * H;
phi = eig(S);
K = 3

[Ptx_K,state] = maxpower_Kstreams(phi,K,'uniform_rate')