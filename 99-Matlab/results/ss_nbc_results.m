clc; clear all; close all;
format long;

% fem
fem1 = load("./NonlinearFEM2D_T6_SS_NBC_2x2.mat");
fem2 = load("./NonlinearFEM2D_T6_SS_NBC_4x4.mat");
fem3 = load("./NonlinearFEM2D_T6_SS_NBC_6x6.mat");
fem4 = load("./NonlinearFEM2D_T6_SS_NBC_8x8.mat");
fem5 = load("./NonlinearFEM2D_T6_SS_NBC_10x10.mat");
fem_dof = [2*size(fem1.model.Nodes,1), 2*size(fem2.model.Nodes,1),...
    2*size(fem3.model.Nodes,1), 2*size(fem4.model.Nodes,1), 2*size(fem5.model.Nodes,1)];
fem_relerrdisp = [fem1.result.Relerrdisp, fem2.result.Relerrdisp,...
    fem3.result.Relerrdisp, fem4.result.Relerrdisp, fem5.result.Relerrdisp];

% csfem
% fem
csfem1 = load("./NonlinearSFEM2D_cell_T6_SS_NBC_2x2.mat");
csfem2 = load("./NonlinearSFEM2D_cell_T6_SS_NBC_4x4.mat");
csfem3 = load("./NonlinearSFEM2D_cell_T6_SS_NBC_6x6.mat");
csfem4 = load("./NonlinearSFEM2D_cell_T6_SS_NBC_8x8.mat");
csfem5 = load("./NonlinearSFEM2D_cell_T6_SS_NBC_10x10.mat");
csfem_relerrdisp = [csfem1.result.Relerrdisp, csfem2.result.Relerrdisp,...
    csfem3.result.Relerrdisp, csfem4.result.Relerrdisp, csfem5.result.Relerrdisp];
