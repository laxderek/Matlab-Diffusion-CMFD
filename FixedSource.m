%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
%clc
close all

solver = Solver();

%%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
solver.MethodGauss = 1;
solver.MethodEigs = 0;
solver.max_iters = 10000;
solver.convergence = 10^-12;
solver.max_iters_coarse = 1000;
solver.convergence_coarse = 10^-12;
solver.gridReductionFactor = [5 1 1];
solver.CMFD = 0;
solver.iterationsBetweenCMFD = 1;
solver.outerprintlevel = 0;
solver.innerprintlevel = 0;
solver.verify = 0;
solver.verifyCoarse = 0;
solver.figures = 1;
solver.fixedSource = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
solver.ng = 2;

%%% First Material
mat = Material();
%D = [1.5];
D = [1.5 0.3];
%sigA = [0.15];
sigA = [0.05 .5];
%nusigF = [0.15];
chi = [1.0 0.0];
nusigF = [0.02 .48];
%sigS = [0.05];
sigS = [[0.05 0.1]
        [0.0 0.5]]; 
mat.D = D;
%tot = [.2];
tot = [.2 1];
mat.chi = chi;
mat.nusigF = nusigF;
mat.sigT = tot;
mat.sigA = sigA;
mat.sigS = sigS;
mat.nusigF = nusigF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%Geometry Definitions%%%%%%%%%%%%%%%%%%
%Dimensions of problem. This value superceeds all other 
%settings and will override them where necessary.
solver.dimensions = 1;
%Size of each dimension in cm
solver.dim = [50 1 1];
%Number of cells in each dimension
x_len = 10; y_len = 1; z_len = 1;
%Set boundary conditions via albedo. For reference:
% 0 = Vacuum
% 1 = Reflective
solver.BC = [1 1 1];
mesh = Mesh(x_len,y_len,z_len,solver.ng);
mesh = mesh.setNumMats(1);
mesh.mats = [mat];
mesh = mesh.setAllMat(1,mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solver.mesh = mesh;
solver = solver.solve();
%solver = solver.calculateCurrents();
%solver.checkBalance();
%close all
 