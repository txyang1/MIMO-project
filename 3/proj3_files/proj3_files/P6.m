clear all
close all
clc

load('exampleBlockDiagMIMOBC.mat')

Ptx = 10;
P = 10.^(Ptx/10);

[ Rsum ] = BlockDiagBCEqualPower( H,C,P )