clc
close all
clear all

MODEL_DIR='../../../..';
addpath([MODEL_DIR '/matlab'])

%% Define output file name
meshID=0; % creates ../N/N_1.txt and ../T/T_1.txt
useRegions=1;
targetElements=1e4;
filename='bicrystal'; % this creates file filename.poly

R=1000;       % radius of prism
H=4*R;        % height of prism
np=8;         % number of points along circumference

%V=pi*R^2*H;    % volume
%N=[1 0 3]';    % normal to GB plane
N=[1 0 1]';
P0=[0 0 100]'; % a point on the GB plane
A=[0 1 1;1 0 1;1 1 0]/sqrt(2); % matrix of primitive lattice vectors for FCC

writeBicrystalPoly(meshID,useRegions,targetElements,filename,R,H,np,N,P0,A)

