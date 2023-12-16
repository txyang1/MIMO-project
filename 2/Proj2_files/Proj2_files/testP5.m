clear all
close all
clc

load('exampleMAC.mat');
w = [1;1];
[Q, Cwsr] = maxWSRmac(H,P,w);