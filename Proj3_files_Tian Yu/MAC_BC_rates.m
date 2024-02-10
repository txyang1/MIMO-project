function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
% function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
%
% The function calculates the achievable rates of the MAC and
% BC for given channels Hk and transmit covariance matrices Qk and Sk and
% a given decoding and encoding order.
%
% Input Specifications:
% H: Kx1 cell array with user channels Hk
% Q: Kx1 cell array with MxM matrix Qk per cell
% S: Kx1 cell array with NxN matrix Sk per cell
% order: 1xK vector in {1,...,K} with the BC encoding order
%
% Output Specifications:
% R_BC:  Kx1 array with achievable rates of the BC
% R_MAC: Kx1 array with achievable rates of the MAC

[M,N] = size(H{1});
K = length(H);

R_BC = zeros(K,1);
R_MAC = zeros(K,1);

I_M = eye(M);
I_N = eye(N);

for k = 1:K
    other_ind = find(order>k);
    BC = I_M;
    if ~isempty(other_ind)
        for i = other_ind
            BC = BC + H{k}*S{i}*H{k}';
        end
    end
    R_BC(k) = log2(det(I_M + BC^(-1)*H{k}*S{k}*H{k}'));
end

for k = K:-1:1
    other_ind = find(order<k);
    MAC = I_N;
    if ~isempty(other_ind)
        for i = other_ind
            MAC = MAC + H{i}'*Q{i}*H{i};
        end
    end
    R_MAC(k) = log2(det(I_N + MAC^(-1)*H{k}'*Q{k}*H{k}));
end

end







