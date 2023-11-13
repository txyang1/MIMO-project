function [ Ptx ] = activeStreams_waterfilling( phi )
% function [ Ptx ] = activeStreams_waterfilling( phi )
%
% Function computes the power values for which the waterfilling power
% allocation switches from K to K+1 active streams
%
% Input:
% phi: vector of eigenmode coefficients phi1,...,phiN
% Output:
% Ptx: vector of power values (size N-1x1) for which the power allocation
% switches from K to K+1 active streams

% TODO

phi = phi(:);
phi_frac = 1./phi;
phi_frac = sort(phi_frac);
Ptx = zeros(length(phi)-1,1);

for K = 1:length(phi)-1
    Ptx(K) = K * phi_frac(K+1) - sum(phi_frac(1:K));
end


end

