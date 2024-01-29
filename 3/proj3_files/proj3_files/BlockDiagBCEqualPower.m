function [ Rsum ] = BlockDiagBCEqualPower( H,C,Ptx )
% function [ Rsum ] = BlockDiagBCEqualPower( H,C,Ptx )
%
% The function computes the maximum achievable rate of a MIMO BC with block
% diagonalization and equal power allocation
%
% Inputs
% H: K x 1 cell array with channel Hk per cell
% C: K x 1 cell array with noise covariance C_nk per cell
% Ptx: joint available transmit power
% Outputs
% Rsum: achievable sum rate

[M,N] = size(H{1});
K = length(H);
R = zeros(K,1);

for k = 1:K
    H_ = [];
    for j = 1:K
        if j~=k
            H_ = [H_;H{j}];
        else
            continue
        end
    end
%     H_
    Pk = eye(N) - H_'*(H_*H_')^(-1)*H_;
    [Nk,Dk] = eig(Pk);
    dk = round(real(diag(Dk)));
%     dk == 1
    Vk = Nk(:,dk==1);
%     Vk = Nk*diag(sqrt(diag(Dk)));
    Xk = Vk'*H{k}'*(C{k})^(-1)*H{k}*Vk;
    [~,Phik] = eig(Xk);
    phi = real(diag(Phik));
    psi = waterfilling(phi,Ptx/K);
%      log2(1+phi.*psi)
    R(k) = real(sum(log2(1+phi.*psi),'all'));
end

Rsum = sum(R,'all');
