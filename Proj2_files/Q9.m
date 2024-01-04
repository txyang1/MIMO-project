clear all
close all
clc

load('exampleMAC.mat');
M = size(H,1);
%% Q4

Q4 = zeros(size(Q));

for i = 1:size(H,3)
    X = H(:,:,i)' * H(:,:,i);
    [Q4(:,:,i), ~] = ratemaxQk(X, P(i));
end

[R4, R4sum] = ratesMAC(Q4, H);

%% Q6

% user1 single-user capacity
Q6_1 = zeros(size(Q));
k = 1;
j = 2;

Xk = H(:,:,k)' * H(:,:,k);
[Q6_1(:,:,k), ~] = ratemaxQk(Xk, P(k));
Xj = H(:,:,j)' * (eye(M) + H(:,:,k)*Q6_1(:,:,k)*H(:,:,k)')^-1 * H(:,:,j);
[Q6_1(:,:,j), ~] = ratemaxQk(Xj, P(j));
[R6_1, R6_1sum] = ratesMAC(Q6_1, H);

% user2 single-user capacity
Q6_2 = zeros(size(Q));
k = 2;
j = 1;

Xk = H(:,:,k)' * H(:,:,k);
[Q6_2(:,:,k), ~] = ratemaxQk(Xk, P(k));
Xj = H(:,:,j)' * (eye(M) + H(:,:,k)*Q6_1(:,:,k)*H(:,:,k)')^-1 * H(:,:,j);
[Q6_2(:,:,j), ~] = ratemaxQk(Xj, P(j));
[R6_2, R6_2sum] = ratesMAC(Q6_2, H);

%% Q8

[Q8, ~, ~] = iterWaterfill(H, P);
[R8, R8sum] = ratesMAC(Q8, H);

%% plot
singleAndSumBounds = figure(1);
hold on
singleAndSumBounds = plotRegionMAC(R4,singleAndSumBounds);
singleAndSumBounds = plotRegionMAC(R6_1,singleAndSumBounds);
singleAndSumBounds = plotRegionMAC(R6_2,singleAndSumBounds);
singleAndSumBounds = plotRegionMAC(R8,singleAndSumBounds);

legend('single user case: both','single user case: 1', 'single user case: 2', 'iterative waterfilling','Location','southwest')
xlabel('R1')
ylabel('R2')





