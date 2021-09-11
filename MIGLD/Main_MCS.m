% This script aims to obtain the extreme valuse of a linear structure under
% non-stationary stochastic seismic ground motions by MCS

clear all; close all;clc

%% MCS 

tic

num = 1e5;
[ gg ] = Non_stationary_seismic_motoins_mcs( num );
[ D1_0] = Newmark_linear_shear_frame_structure( num, gg );

toc


figure(1);
[jjj,iii]=hist(D1_0,60); 
jjj=jjj/length(D1_0)/mean(diff(iii));
b=bar(iii,jjj,1);