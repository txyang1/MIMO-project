clear all;
clc
load('example_channels.mat');
S = H' * inv(Cn) * H;
phi = eig(S);

Ptx = activeStreams_waterfilling(phi);

[psi,mu,K] = waterfilling(phi,100);