clear all
close all
clc

load('exampleMIMOBC.mat');
Ptx = -15:5:30;
[fig] = plotSumRateBC(H,Ptx);