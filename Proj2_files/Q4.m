clear all
close all
clc

load('exampleMAC.mat');
for i = 1:size(H,3)
    X = H(:,:,i)' * H(:,:,i);
    [Q, C] = ratemaxQk(X, P(i))
end