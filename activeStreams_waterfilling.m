function [ Ptx ] = activeStreams_waterfilling( phi )
%
% Function computes the power values for which the waterfilling power
% allocation switches from K to K+1 active streams
%
% Input:
% phi: vector of eigenmode coefficients phi1,...,phiN
% Output:
% Ptx: vector of power values (size N-1x1) for which the power allocation
% switches from K to K+1 active streams

phi = sort(phi,'descend'); %sort it in non-increasing manner
inv_phi = 1./phi; %invert the vector phi 
N = numel(phi); %get number of elements in vector phi
K = 1:(N-1);
Ptx = zeros(N-1,1); %ready to save data

for i=1:(N-1)
    Ptx(i) = K(i)*inv_phi(K(i)+1)-sum(inv_phi(1:i)); %loop
end

end

%Team members: Tian Yu, Tingxin Yang