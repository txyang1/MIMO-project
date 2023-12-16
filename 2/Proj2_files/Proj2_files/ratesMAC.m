function [R,Rsum] = ratesMAC(Q,H)
%
% The function ratesMAC computes the rate region coordinates of
% the MIMO MAC for given transmit covariance matrices.
%
% Inputs
% H: M x N x 2 array of the given channel realizations, where
%    N is the number transmit antennas
%    M is the number receive antennas
% Q: M x N x 2 array of the given tranmit covariance matrices realizations
% Outputs
% R: 2 x 2 matrix of rate region coordinates R = [r_A, r_B]
% Rsum: achievable sum-rate

M = size(H,1);

% TODO
Q1 = Q(:,:,1);
Q2 = Q(:,:,2);
H1 = H(:,:,1);
H2 = H(:,:,2);

Rsum = real(log2(det(eye(M) + H1*Q1*H1' + H2*Q2*H2')));
C1 = real(log2(det(eye(M) + H1*Q1*H1')));
C2 = real(log2(det(eye(M) + H2*Q2*H2')));
rA = [Rsum-C2, C2].';
rB = [C1, Rsum-C1].';
R = [rA, rB];
