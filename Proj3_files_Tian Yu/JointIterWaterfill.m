function [Q,Csum,Rsum] = JointIterWaterfill(H,Ptx,epsilon)
% Function [Q,Csum,Rsum] = JointIterWaterfill(H,E,epsilon)
% 
% The function calculates the sum capacity and the corresponding
% transmit covariance matrices of the K user MIMO MAC with a sum
% transmit power constraint via a joint iterative waterfilling.
%
% Inputs
% H: M x N x 2 array of users' channel matrices
% Etx: joint available transmit power
% epsilon: accuracy of the algorithm (sopping criteria)
% Outputs
% Q: N x N x 2 array of users' optimal transmit covariance matrices
% Csum: sum capacity of the MIMO MAC
% Rsum: vecotr of sum rate values for the iterations


if nargin<3, epsilon = 1e-6; end
maxiter = 1e4;

K = length(H);
[M,N] = size(H{1});


H1 = H{1};
H2 = H{2};

Q1 = Ptx/(2*M)*eye(M);
Q2 = Ptx/(2*M)*eye(M);

Q = cell(K,1);

ii=0;

while true
    ii=ii+1;
    
    X1 = (H1/(H2'*Q2*H2 + eye(N)))*H1';
    X2 = (H2/(H1'*Q1*H1 + eye(N)))*H2';
    
    [V1,Phi1] = eig(X1);
    [V2,Phi2] = eig(X2);
    
    phi1 = diag(Phi1);
    phi2 = diag(Phi2);
    
    psi = waterfilling(abs([phi1(:);phi2(:)]), Ptx);
    
    psi1 = psi(1:M);
    psi2 = psi(M+1:end);
    
    Q1_tmp = V1*diag(psi1)*V1';
    Q2_tmp = V2*diag(psi2)*V2';
    %Q1_tmp = 0.5*(Q1_tmp+Q1_tmp');
    %Q2_tmp = 0.5*(Q2_tmp+Q2_tmp');
    
    Rsum(ii) = real(log2(det( eye(N) + H1'*Q1_tmp*H1 + H2'*Q2_tmp*H2)));
    %norm(Q1_tmp-Q1,'fro')^2 + norm(Q2_tmp-Q2,'fro')^2
    if(norm(Q1_tmp-Q1,'fro')^2 + norm(Q2_tmp-Q2,'fro')^2 <= epsilon)
        Q{1} = Q1_tmp;
        Q{2} = Q2_tmp;
        break;
    end
    
    Q1 = Q1_tmp;
    Q2 = Q2_tmp;
    
end

Csum = Rsum(end);
