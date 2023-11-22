function [ Ptx ] = activeStreams_waterfilling( phi )
%
% Function computes the power values for which the waterfilling power
% allocation switches from K to K+1 active streams
%
% Input:
phi = sort(phi,'descend');
phi_inv = 1./phi;

% phi: vector of eigenmode coefficients phi1,...,phiN
% Output:
% Ptx: vector of power values (size N-1x1) for which the power allocation

Ptx = zeros(length(phi)-1,1);
for K = 1: length(phi)-1
    Ptx(K) = K*phi_inv(K+1)-sum(phi_inv(1:K));
end
% switches from K to K+1 active streams

% TODO

end





