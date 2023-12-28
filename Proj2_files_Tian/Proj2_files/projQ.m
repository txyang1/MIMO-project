function [Q] = projQ(Qp,P)
% Function [Q] = projQ(Qp,P)
%
% Orthogonal projections of the N x N Hermitian 
% non-negative definite matrices in Qp onto the 
% set of positive semidefinite matrices satisfying trace(Q)<=P.
%
% Inputs
% Qp: N x N x 2 array of the given Hermitian non-negative definite matrices
% P: 2 x 1 vector of transmit powers
% Outputs
% Q: N x N x 2 array of positive semidefinite projected matrices Q1,Q2


N = size(Qp,2);
K = size(Qp,3);
Q = zeros(N,N,K);

%TODO
Q1_dash = Qp(:,:,1);
Q2_dash = Qp(:,:,2);


%This code is also strange. Problems?


if trace(Q1_dash)<=P(1) 
    Q(:,:,1) = Q1_dash;
else
    [U1,phi1] = eig(Q1_dash);
    psi1 = waterspilling(diag(phi1),P(1));
    Q(:,:,1) = norm(Q1_dash-U1*diag(psi1)*U1','fro').^2;
end

if trace(Q2_dash)<=P(2) 
    Q(:,:,2) = Q2_dash;
else
    [U2,phi2] = eig(Q2_dash);
    psi2 = waterspilling(diag(phi2),P(2));
    Q(:,:,2) = norm(Q2_dash-U2*diag(psi2)*U2','fro').^2;
end

%Team members: Tingxin Yang, Tian Yu
