% test_findPositionFib.m
%
%
%
% Jorge Llop - Apr 13, 2018
clear all;
close all;

addpath(genpath('utils'));

hcstt_Initialize(true);
hcstt_NewFlatForDM('ImageSharpening_fminconIt2_Apr1');
[actxc_fib,ang_fib] = hcstt_FindPosiotionFiberv4_coarse(2.4,0.0)
hcstt_DisconnectDevices();
