clear all;
clc
load('example_channels.mat');
S = H' * inv(Cn) * H;
phi = eig(S);

[Ptx] = activeStreams_mmse(phi);

[psi,mu,K] = mmseallocation(phi,25);