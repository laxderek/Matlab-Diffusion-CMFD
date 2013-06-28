%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
clc
close all

solver = Solver();

%%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
solver.MethodGauss = 1;
solver.MethodEigs = 0;
solver.max_iters = 10000;
solver.convergence = 10^-12;
solver.max_iters_coarse = 1000;
solver.convergence_coarse = 10^-12;
solver.gridReductionFactor = [2 1 1];
solver.CMFD = 1;
solver.iterationsBetweenCMFD = 1;
solver.outerprintlevel = 1;
solver.innerprintlevel = 0;
solver.verify = 0;
solver.verifyCoarse = 0;
solver.figures = 0;
solver.fixedSource = 0;
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

%%% Second Material
mat2 = Material();
D = [1.0 .1];
sigA = [0.5 .5];
nusigF = [0 0.0];
chi = [0.0 0.0];
sigS = [[0.5 0.1]
        [0.01 4.5]]; 
tot = [1.1 5.01];
mat2.D = D;
mat2.nusigF = nusigF;
mat2.chi = chi;
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
solver.dim = [200 1 1];
%Number of cells in each dimension
x_len = 100; y_len = 1; z_len = 1;
%Set boundary conditions via albedo. For reference:
% 0 = Vacuum
% 1 = Reflective
solver.BC = [0 1 1];
mesh = Mesh(x_len,y_len,z_len,solver.ng);
mesh = mesh.setNumMats(1);
mesh.mats = [mat mat2];
mesh = mesh.setAllMat(1,mat);
%mesh = mesh.setMatAtLoc(1,1,1,2,mat2);
%mesh = mesh.setMatAtLoc(2,1,1,2,mat2);
%mesh = mesh.setMatAtLoc(3,1,1,2,mat2);
%mesh = mesh.setMatAtLoc(4,1,1,2,mat2);
%mesh = mesh.setMatAtLoc(5,1,1,2,mat2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

solver.mesh = mesh;
solver = solver.solve();
k1 = solver.k;
%solver = solver.calculateCurrents();
%solver.checkBalance();
%close all

solver.CMFD = 0;
solver = solver.solve();
k2 = solver.k;

disp(k1)
disp(k2)
disp(k1-k2)