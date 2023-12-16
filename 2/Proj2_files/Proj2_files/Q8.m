clear all
close all
clc

load('exampleMAC.mat')

[Q, Csum, Rsum] = iterWaterfill(H, P);