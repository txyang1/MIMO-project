clear all
close all
clc

load('exampleMIMOBC.mat');
% H1 = cell(size(H));
% for i = 1:length(H1)
%     H1{i} = rand(size(H{i}));
% end

Ptx = -15;
P = 10.^(Ptx/10);

[ Q, Csum ] = DualMACSumRateMaximization( H, P );

