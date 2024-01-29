clear all
close all
clc

load('exampleMIMOBC.mat');
order = [1,2]';
S{1} = eye(3);
S{2} = 2*eye(3);

[R_BC, R_MAC] = MAC_BC_rates(H,Q,S,order);