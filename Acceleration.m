%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
clc
close all

solver = Solver();

%%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
solver.MethodGauss = 1;
solver.MethodEigs = 0;
solver.max_iters = 1000;
solver.convergence = 10^-11;
solver.max_iters_coarse = 100;
solver.convergence_coarse = 10^-6;
solver.gridReductionFactor = [2 1 1];
solver.CMFD = 1;
solver.iterationsBetweenCMFD = 1;
solver.outerprintlevel = 1;
solver.innerprintlevel = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
solver.ng = 2;

%%% First Material
mat = Material();
D = [1.5 0.3];
sigA = [0.05 .5];
nusigF = [[0.02 0.0]
          [0.48 0.0]];
sigS = [[0.05 0.1]
        [0.0 0.5]]; 
mat.D = D;
tot = [.2 1];
mat.nusigF = nusigF;
mat.sigT = tot;
mat.sigA = sigA;
mat.sigS = sigS;
mat.nusigF = nusigF;

%%% Second Material
mat2 = Material();
D = [1.0 .1];
sigA = [0.5 0];
nusigF = [[0 0.0]
          [0 0.0]];
sigS = [[0.5 0.1]
        [0.01 5]]; 
tot = [1.1 5.01];
mat2.D = D;
mat2.nusigF = nusigF;
mat2.sigT = tot;
mat2.sigA = sigA;
mat2.sigS = sigS;
mat2.nusigF = nusigF;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Geometry Definitions%%%%%%%%%%%%%%%%%%
%Dimensions of problem. This value superceeds all other 
%settings and will override them where necessary.
solver.dimensions = 1;
%Size of each dimension in cm
solver.dim = [10 1 1];
%Number of cells in each dimension
x_len = 6; y_len = 1; z_len = 1;
%Set boundary conditions via albedo. For reference:
% 0 = Vacuum
% 1 = Reflective
solver.BC = [0 1 1];
mesh = Mesh(x_len,y_len,z_len,solver.ng);
mesh = mesh.setNumMats(1);
mesh.mats = [mat];
mesh = mesh.setAllMat(1,mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solver.mesh = mesh;
solver.solve();
 
%close all
 