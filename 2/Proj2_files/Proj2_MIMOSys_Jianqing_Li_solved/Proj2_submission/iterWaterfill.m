function [Q,Csum,Rsum] = iterWaterfill(H,P,epsilon)
% 
% The function calculates the sum capacity and the corresponding
% transmit covariance matrices of the K user MIMO MAC via
% the iterative waterfilling principle.
%
% Inputs
% H: M x N x 2 array of users' channel matrices
% P: 2 x 1 vector of users' available transmit powers
% epsilon: accuracy of the algorithm (sopping criteria)
% Outputs
% Q: N x N x 2 array of users' optimal transmit covariance matrices
% Csum: sum capacity of the MIMO MAC
% Rsum: vector of sum rate values for the iterations
%

if nargin<3, epsilon = 1e-4; end
maxiter = 1e4;
M = size(H,1);
N = size(H,2);

% ToDo: Initialize Q1 and Q2
Q1 = zeros(N);
Q2 = zeros(N);

ii=0;
while(ii <= maxiter)
    ii=ii+1;
    
    % ToDo: Determine update for Q1 and Q2 and achievable sum rate after
    % the update
    
    [Q1_tmp, ~] = ratemaxQk(H(:,:,1)'*(eye(M)+H(:,:,2)*Q2*H(:,:,2)')^-1*H(:,:,1), P(1));
    [Q2_tmp, ~] = ratemaxQk(H(:,:,2)'*(eye(M)+H(:,:,1)*Q1*H(:,:,1)')^-1*H(:,:,2), P(2));
    
    Rsum(ii) = real(log2(det(eye(M) + H(:,:,1)*Q1_tmp*H(:,:,1)' + H(:,:,2)*Q2_tmp*H(:,:,2)')));
    

    %ToDo: Insert converegence criterion
    if(norm(Q1_tmp-Q1, 'fro')+norm(Q2_tmp-Q2, 'fro') < epsilon)
        break;
    end
    
    % Update Q1 and Q2 for next iteration
    Q1 = Q1_tmp;
    Q2 = Q2_tmp;
end

Q(:,:,1) = Q1_tmp;
Q(:,:,2) = Q2_tmp;
Csum = Rsum(end);
