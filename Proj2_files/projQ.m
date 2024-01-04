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
Qp1 = Qp(:,:,1);
Qp2 = Qp(:,:,2);

phi1 = real(eig(Qp1));
[U1,~] = eig(Qp1);
[psi1,~,~] = waterspilling(phi1,P);
Q(:,:,1) = U1*diag(psi1)*U1';

phi2 = real(eig(Qp2));
[U2,~] = eig(Qp2);
[psi2,~,~] = waterspilling(phi2,P);
Q(:,:,2)= U2*diag(psi2)*U2';

%Team members: Tingxin Yang, Tian Yu



