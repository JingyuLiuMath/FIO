clear;
close all;

data_path = "/scratch/jyliu/FIO/2d_const_amplitude/data/";

originalPath = path;
addpath('../../extern/FastBF.m/src');

exp_phi_func = @(x, xi) fun_2D(x, xi);

p_list = (6 : 9)';
num_n = length(p_list);

r_bf = 10;
tol_bf = 1e-8;

min_points = 256;
tol_hss_list = [1e-3];
num_tol_hss = length(tol_hss_list);

tol_cg = 1e-12;
maxit_cg = 500;

num_sample = 256;