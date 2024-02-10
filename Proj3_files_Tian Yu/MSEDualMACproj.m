function [Qblk] = MSEDualMACproj(QblkHat, Ptx)
% function [Qblk] = MSEDualMACproj(QblkHat, Ptx)
%
% MSEDualMACproj computes the orthogonal projection of a block diagonal
% matrix QlkbHat such that that Qrpoj is block diagonal, positive
% semi-definite and fulfilles the power consraint with the sum power Ptx.
% Input:
% QblkHat: KM x KM block diagonal matrix with transmit covariance matrices
%          in each block
% Ptx: Available sum transmit power
% Output:
% Qblk: KM x KM block diagonal matrix that is the orthogonal projection of
%       QblkHat onto the constraint set

QHat = mat2cell(QblkHat,[2 2],[2 2]);

[U1,phi1] = eig(QHat{1,1});
psi1 = waterspilling(phi1,Ptx);
Psi1 = diag(psi1);
Q1 = U1*Psi1*U1';

[U2,phi2] = eig(QHat{2,2});
psi2 = waterspilling(phi2,Ptx);
Psi2 = diag(psi2);
Q2 = U2*Psi2*U2';

Q = cell(2,1);
Q{1} = Q1;
Q{2} = Q2;

Qblk = blkdiag(Q{:});
end

