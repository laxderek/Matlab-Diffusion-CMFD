%Derek Lax
%1D/2D/3D Finite Volume Solver with CMFD


clear
clc
close all

%%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
MethodGauss = 1;
MethodEigs = 0;
max_iters = 10000;
convergence = 10^-10;
gridReductionFactor = [4 1 1];
iterationsBetweenCMFD = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
ng = 2;
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
dimensions = 1;
%Size of each dimension in cm
dim = [2000 1 1];
%Number of cells in each dimension
x_len = 100; y_len = 1; z_len = 1;
%Set boundary conditions via albedo. For reference:
% 0 = Vacuum
% 1 = Reflective
BC = [1 1 1];
mesh = Mesh(x_len,y_len,z_len,ng);
mesh = mesh.setNumMats(1);
mesh.mats = [mat];
mesh = mesh.setAllMat(1,mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   END USER INPUT                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%Build Problem Geometry%%%%%%%%%%%%%%%%%%%%
%Evenly space mesh points
if (dimensions > 3 || dimensions < 1)
    disp('Invalid number of dimensions');
    break
end

% if (dimensions < 3)
%    if (z_len ~= 1)
%        disp('Invalid number of mesh cells in z-axis.');
%        disp('Less than 3D. Number of cells in z-axis must be set to 1');
%        break
%    end
%    dim(3) = 1;
%    BC(3) = 1;
% end
% 
% if (dimensions < 2)
%    if (y_len ~= 1)
%        disp('Invalid number of mesh cells in y-axis.');
%        disp('Problem is 1D. Number of cells in y-axis must be set to 1');
%        break
%    end
%    dim(2) = 1;
%    BC(2) = 1;
% end

mesh = mesh.setEqualDist(1,dim(1),dim(2),dim(3));
%Set boundary conditions
mesh = mesh.setBC(BC(1), BC(2), BC(3));%set BC
%Compute total number of mesh cells
total_mesh = mesh.x * mesh.y * mesh.z * mesh.g;


%%%%%%%%%%%%%%%%%%%%%%Run Solver%%%%%%%%%%%%%%%%%%%%%%%%%%

[M mesh] = BuildLossMatrix(mesh, 1, false, 0);
F = BuildProductionMatrix(mesh, 1);

%Initial Guess:
format long
k = 1;
kold = .1;
phiold = ones(total_mesh,1);
if MethodGauss == 1
    for iters = 1:max_iters

        if (abs(k-kold) < convergence)
            break;
        end
        kold = k;
        b = F*phiold/kold;
        phi = M\b;
        k = kold * sum(phi)/sum(phiold);
         phi = phi / sum(sum(phi));
        phiold = phi;

    end
elseif MethodEigs == 1
    k = eigs(M\F);
end



ktrue = k;
disp k;
disp (k);
disp('iters');
disp(iters);


INF = [1 mesh.g mesh.x mesh.y mesh.z];
%calculate all surface currents
for i1 = 1:mesh.x
    for j1 = 1:mesh.y
        for k1 = 1:mesh.z
            for g1 = 1:ng
                for surf = 1:6
                    %check for boundary
                    xyz = round((surf / 2)); %1=x;2=y;3=z
                    %%determine left or right surface
                    side = -2 * mod(surf,2) + 1; % -1=left, 1=right
                    %%determine neighbor index
                    %if doesn't exist (boundary), returns -1
                    %if doesn't exist (boundary), returns -1
                    if (xyz == 1)
                        i2 = i1 + side;
                        j2 = j1;
                        k2 = k1;
                    elseif (xyz == 2)
                        i2 = i1;
                        j2 = j1 + side;
                        k2 = k1;
                    else %(xyz == 3) 
                        i2 = i1;
                        j2 = j1;
                        k2 = k1 + side;		        	
                    end
                    neigh = checkBoundary(i1,j1,k1,mesh.x,mesh.y,mesh.z,xyz,side);
                    if (neigh == 1)
                    %not boundary
                        %calculate current for each boundary
                        mesh.Jsurf(i1,j1,k1,g1,surf) = -mesh.dtilde(i1,j1,k1,g1,surf) * ...
                            (phi(indexToMat(i1,j1,k1,g1,INF)) - phi(indexToMat(i2,j2,k2,g1,INF)));                                
                    else
                    %boundary
                    %Only works for vacuum currently
                        mesh.Jsurf(i1,j1,k1,g1,surf) = -mesh.dtilde(i1,j1,k1,g1,surf) * ...
                            phi(indexToMat(i1,j1,k1,g1,INF));
                    end
                    if (side == 1)
                            mesh.Jsurf(i1,j1,k1,g1,surf) = -mesh.Jsurf(i1,j1,k1,g1,surf);
                    end
                end
            end
        end
    end
end    

%save flux
mesh.phi = phi;
%Make coarse mesh from fine mesh
coarse = Mesh.fineToCoarse(mesh,phi,gridReductionFactor,dimensions);

[M2 coarse] = BuildLossMatrix(coarse, 1, true, mesh, 1.0);
F2 = BuildProductionMatrix(coarse, 1);
 
total_mesh = coarse.x * coarse.y * coarse.z * coarse.g;
%Initial Guess:
format long
k2 = 1;
kplot = zeros(max_iters,1);
kold2 = .1;
phiold = ones(total_mesh,1);
if MethodGauss == 1
    for iters = 1:max_iters

        if (abs(k2-kold2) < convergence)
            break;
        end
        kold2 = k2;
        b2 = F2*phiold/kold2;
        phi2 = M2\b2;
        k2 = kold2 * sum(phi2)/sum(phiold);
        kplot(iters) = k2;
        phi2 = phi2 / sum(sum(phi2));
        phiold = phi2;
        flux2 = reshape(phi2,ng,numel(phi2)/ng,1,1);
        
     end
 elseif MethodEigs == 1
     k2 = eigs(M2\F2);
end
disp k;
disp (k2);
disp('iters');
disp(iters);
 
close all



% 
% solver = Solver();
% 
% %%%%%%%%%%%%Problem Solving Specifications%%%%%%%%%%%%%%%%
% solver.MethodGauss = 1;
% solver.MethodEigs = 0;
% solver.max_iters = 10000;
% solver.convergence = 10^-10;
% solver.max_iters_coarse = 10000;
% solver.convergence_coarse = 10^-12;
% solver.gridReductionFactor = [2 1 1];
% solver.CMFD = 0;
% solver.iterationsBetweenCMFD = 1;
% solver.outerprintlevel = 0;
% solver.innerprintlevel = 0;
% solver.verify = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%% Material Defintions %%%%%%%%%%%%%%%%%%%%%%
% solver.ng = 2;
% 
% mat = Material();
% %D = [1.5];
% D = [1.5 0.3];
% %sigA = [0.15];
% sigA = [0.05 .5];
% %nusigF = [0.15];
% chi = [1.0 0.0];
% nusigF = [0.02 .48];
% %sigS = [0.05];
% sigS = [[0.05 0.1]
%         [0.0 0.5]]; 
% mat.D = D;
% %tot = [.2];
% tot = [.2 1];
% mat.chi = chi;
% mat.nusigF = nusigF;
% mat.sigT = tot;
% mat.sigA = sigA;
% mat.sigS = sigS;
% mat.nusigF = nusigF;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%Geometry Definitions%%%%%%%%%%%%%%%%%%
% %Dimensions of problem. This value superceeds all other 
% %settings and will override them where necessary.
% solver.dimensions = 1;
% %Size of each dimension in cm
% solver.dim = [200 1 1];
% %Number of cells in each dimension
% x_len = 10; y_len = 1; z_len = 1;
% %Set boundary conditions via albedo. For reference:
% % 0 = Vacuum
% % 1 = Reflective
% solver.BC = [0 1 1];
% mesh = Mesh(x_len,y_len,z_len,solver.ng);
% mesh = mesh.setNumMats(1);
% mesh.mats = [mat];
% mesh = mesh.setAllMat(1,mat);
% %mesh.setMatAtLoc(self,1,1,1,mat2,material) 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver.figures = 0;
% solver.mesh = mesh;
% 
% solver = solver.solve();
%
%solver = solver.accelerate();

