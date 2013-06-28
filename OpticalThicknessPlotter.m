%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
clc
close all

max = 50;
min = 1;
step = 10;
num = (max - min) / step;
sigMult = linspace(min,max,num);
opticalThickness = zeros(num,1);
iters = zeros(num,1);
k = zeros(num,1);
iters2 = zeros(num,1);
k2 = zeros(num,1);
for i = 1:num
    disp(i);
    solver = Solver();

    %%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
    solver.MethodGauss = 1;
    solver.MethodEigs = 0;
    solver.max_iters = 10000;
    solver.convergence = 10^-8;
    solver.max_iters_coarse = 1000;
    solver.convergence_coarse = 10^-8;
    solver.gridReductionFactor = [2 1 1];
    solver.CMFD = 1;
    solver.iterationsBetweenCMFD = 1;
    solver.outerprintlevel = 1;
    solver.innerprintlevel = 0;
    solver.verify = 0;
    solver.figures = 0;
solver.fixedSource = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
    solver.ng = 2;

    %%% First Material
    mat = Material();
    sigA = [0.05 .5*sigMult(i)];
    nusigF = [0.02 0.48];
    chi = [1.0 0.0];
    sigS = [[0.05 0.1]
            [0.0 0.5]]; 
    tot = [.2 1*sigMult(i)];
    D = [1.5 0.3];
    mat.D = D;
    mat.nusigF = nusigF;
    mat.chi = chi;
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
    solver.dim = [200 1 1];
    %Number of cells in each dimension
    x_len = 100; y_len = 1; z_len = 1;
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
    iters(i) = solver.iters;
    k(i) = solver.k;
    solver.CMFD = 0;
    solver.solve();
    iters2(i) = solver.iters();
    k2(i) = solver.k;
    opticalThickness(i) = 2*tot(2)*solver.dim(1)/x_len;
end 


figure;
plot(opticalThickness,iters);
axis([0 20 0 15]);
xlabel('Optical Thickness (#mfp)');
ylabel('Number of Iterations');