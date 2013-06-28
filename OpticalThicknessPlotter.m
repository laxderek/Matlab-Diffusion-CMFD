%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
clc
close all

for sigMult = 1:10:.25
    
    solver = Solver();

    %%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
    solver.MethodGauss = 1;
    solver.MethodEigs = 0;
    solver.max_iters = 10000;
    solver.convergence = 10^-14;
    solver.max_iters_coarse = 1000;
    solver.convergence_coarse = 10^-14;
    solver.gridReductionFactor = [10 1 1];
    solver.CMFD = 1;
    solver.iterationsBetweenCMFD = 1;
    solver.outerprintlevel = 0;
    solver.innerprintlevel = 0;
    solver.verify = 0;
    solver.figures = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
    solver.ng = 2;

    %%% First Material
    mat = Material();
    sigA = [0.05 .5];
    nusigF = [[0.02 0.0]
              [0.48 0.0]];
    sigS = [[0.05 0.1]
            [0.0 0.5]]; 
    tot = [.2 1];
    D = 1 / 3 * tot;
    mat.D = D;
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
    solver.dim = [100 1 1];
    %Number of cells in each dimension
    x_len = 20; y_len = 1; z_len = 1;
    %Set boundary conditions via albedo. For reference:
    % 0 = Vacuum
    % 1 = Reflective
    solver.BC = [0 1 1];
    mesh = Mesh(x_len,y_len,z_len,solver.ng);
    mesh = mesh.setNumMats(1);
    mesh.mats = [mat];
    mesh = mesh.setAllMat(1,mat);
    %mesh.setMatAtLoc(self,1,1,1,mat2,material) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    solver.mesh = mesh;
    solver = solver.solve();
end 