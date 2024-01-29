clear all
close all
clc

load('exampleMIMOBC.mat');
Ptx = 10;
[ Q, Csum ] = DualMACSumRateMaximization( H, Ptx );

order = [1,2]
[S] = MACtoBCtransform(Q,H,order);

S1 = S{1}
trace(S1)

S2 = S{2}
trace(S2)