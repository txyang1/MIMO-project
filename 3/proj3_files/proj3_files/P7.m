clear all
close all
clc

load('exampleBlockDiagMIMOBC.mat')

Ptx = 40;
P = 10.^(Ptx/10);

[ Rsum ] = BlockDiagBC( H,C,P )