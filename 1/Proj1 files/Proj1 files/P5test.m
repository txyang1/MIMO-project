clear all;
clc
load('example_channels.mat');
S = H' * inv(Cn) * H;
phi = eig(S);
phi
[psi,K] = uniform_rate(phi,48)