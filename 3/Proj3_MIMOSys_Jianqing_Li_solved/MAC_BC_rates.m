function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
% function [ R_BC, R_MAC ] = MAC_BC_rates( H, Q, S, order )
%
% The function calculates the achievable rates of the MAC and
% BC for given channels Hk and transmit covariance matrices Qk and Sk and
% a given decoding and encoding order.
%
% Input Specifications:
% H: Kx1 cell array with channel Hk per cell
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

% S_BC = 0;
% for i = K:-1:1
%     current = find(order==i);
%     pre = find(order==i+1);
%     if ~isempty(pre)
%         S_BC = S_BC + S{pre};
%         R_BC(current) = log2(det(I_N + H{current}'*(I_M+H{current}*S_BC*H{current}')^(-1)*H{current}*S{current}));
%     else
% %        R_BC(current) = log2(det(I_N + (H{current}'*H{current})*S{current}));
%         R_BC(current) = log2(det(I_N + H{current}'*(I_M+H{current}*S_BC*H{current}')^(-1)*H{current}*S{current}));
%     end
% end
% 
% 
% X_MAC = 0;
% for i = 1:K
%     current = find(order==i);
%     pre = find(order==i-1);
%     if ~isempty(pre)
%         X_MAC = X_MAC + H{pre}'*Q{pre}*H{pre};
%         R_MAC(current) = log2(det(I_M + H{current}*(I_N+X_MAC)^(-1)*H{current}'*Q{current}));
%     else
% %        R_MAC(current) = log2(det(I_M + H{current}*H{current}'*Q{current}));
%         R_MAC(current) = log2(det(I_M + H{current}*(I_N+X_MAC)^(-1)*H{current}'*Q{current}));
%     end
%     
% end


for k = 1:K
    current = find(order==k);
    x_index = find(order>k);
    Fk = I_M;
    if ~isempty(x_index)
        for index = x_index
            Fk = Fk + H{current}*S{index}*H{current}';
        end
    end
    R_BC(current) = log2(det(I_M + Fk^(-1)*H{current}*S{current}*H{current}'));
end

for k = K:-1:1
    current = find(order==k);
    x_index = find(order<k);
    Xk = I_N;
    if ~isempty(x_index)
        for index = x_index
            Xk = Xk + H{index}'*Q{index}*H{index};
        end
    end
    R_MAC(current) = log2(det(I_N + Xk^(-1)*H{current}'*Q{current}*H{current}));
end


    
R_BC = real(R_BC);
R_MAC = real(R_MAC);







