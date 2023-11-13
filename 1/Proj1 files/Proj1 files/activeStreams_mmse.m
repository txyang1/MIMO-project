function [ Ptx ] = activeStreams_mmse( phi )
% function [ Ptx ] = activeStreams_mmse( phi )
%
% Function computes the power values for which the MMSE power
% allocation switches from K to K+1 active streams
%
% Input:
% phi: vector of eigenmode coefficients phi1,...,phiN
% Output:
% Ptx: vector of power values (size N-1x1) for which the power allocation
% switches from K to K+1 active streams

% TODO

phi = phi(:);
phi_frac = 1 ./ phi;
phi_frac = sort(phi_frac);
N = length(phi);

Ptx = zeros(N-1,1);
for K = 1:N-1
    Ptx(K) = sqrt(phi_frac(K+1))*sum(sqrt(phi_frac(1:K))) - sum(phi_frac(1:K));
end


end

