clear all

load('example_channels.mat');

S = H' * inv(Cn) * H;
phi = eig(S);

Ptx = activeStreams_waterfilling(phi);

%[psi,mu,K] = waterfilling(phi,Ptx);
%[psi,K] = uniform_rate(phi,Ptx);
%[ psi ] = zf_mmseallocation( phi, Ptx );
%[Ptx_K,state] = maxpower_Kstreams(phi,K,'uniform_rate');
