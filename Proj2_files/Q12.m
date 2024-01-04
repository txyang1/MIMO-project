clear all
close all
clc

load('exampleMAC.mat');
w1 = [1/4,3/4]';
w2 = [3/4,1/4]';


[Cwsum,Q1] = maxWSRmac(H,P(1),w1)

[r1,~] = ratesMAC(Q1,H)
[Cwsum2,Q2] = maxWSRmac(H,P(2),w2)
[r2,~] = ratesMAC(Q2,H)
   
