function [R] = ParetoBound(H,P,S)
% Function [R] = ParetoBound(H,P,S)
%
% The function calculates S coordinate points at the
% Pareto Boundary of the two-user MIMO MAC capacity region
% via varying the weights w_1 and w_2.
%
% Inputs
% H: M x N x 2 array of given channels to the users
% P: 2 x 1 column vector of available transmit power
% S: number of Pareto Boundary sample points S
% Outputs
% R: 2 x 2S matrix of rate region coordinates

M = size(H,1);
N = size(H,2);

R = zeros(2, 2*S);
   
    for s = 1:2:2*S
        w1 = ((s-1)/2)/(2*(S-1));
        w2 = 1 - w1;

        if w1 == 0 || w2 == 0
            
            Q4 = zeros(2,2,2);

            for i = 1:size(H,3)
                X = H(:,:,i)' * H(:,:,i);
                [Q4(:,:,i), ~] = ratemaxQk(X, P(i));
            end
%Team members: Tingxin Yang, Tian Yu
            [R1, ~] = ratesMAC(Q4, H);
            R(:,s:s+1)=R1;
        
        
        elseif w1 == w2 
              [Q8, ~, ~] = iterWaterfill(H, P);
              [R2, ~] = ratesMAC(Q8, H);
              R(:,s:s+1)=R2;
        else
            [~,Q9] = maxWSRmac(H,P,[w1;w2]);
            [R3, ~] = ratesMAC(Q9, H);
            R(:,s:s+1)=R3;
        end
    end

end

       



