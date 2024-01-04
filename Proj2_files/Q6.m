clear all
close all
clc

load('exampleMAC.mat');

k = 1;
j = 2;
N = size(Xk,1);
Xk = H(:,:,k)' * H(:,:,k);
[Qk, Ck] = ratemaxQk(Xk, P(k))
Xj = H(:,:,j)' * (eye(N) + H(:,:,k)*Qk*H(:,:,k)')^-1 * H(:,:,j);
[Qj, Cj] = ratemaxQk(Xj, P(j))