classdef Solver

    properties
	mesh;
    BC;
    dim;
    CMFD;
    M;
    F;
    k;
    ng;
    phi;
    flux;
    iters;
    coarseM;
    coarseF;
    MethodGauss;
    MethodEigs;
    max_iters;
    max_iters_coarse;
    convergence;
    convergence_coarse;
    gridReductionFactor;
    iterationsBetweenCMFD;
    total_mesh;
    dimensions;
    outerprintlevel;
    innerprintlevel;
    end
    
    methods
        function self=Solver() 
            
        end
        function self=solve(self)
          
            
          self = self.checkParams();
          
          
          [self.M self.mesh] = BuildLossMatrix(self.mesh, 1, false, 0);
          self.F = BuildProductionMatrix(self.mesh, 1);
          %Initial Guess:
          format long
          self.k = 1;
          kold = .1;
          phiold = ones(self.total_mesh,1);

          if (self.MethodGauss == 1)
            for i = 1:self.max_iters
              
              if (self.outerprintlevel == 1)
                  text = strcat('Fine Mesh Iteration: ',num2str(i),'. Error in k: ',num2str(abs(self.k-kold)));
                  disp(text);
              end
              if (abs(self.k-kold) < self.convergence)
                  self.iters = i;
                break;
              end
              kold = self.k;
              b = self.F*phiold/kold;
              self.phi = self.M\b;
              self.k = kold * sum(self.phi)/sum(phiold);
              self.phi = self.phi / sum(sum(self.phi));
              if (self.CMFD == 1 && mod(i,self.iterationsBetweenCMFD) == 0) 
                 self.flux = reshape(self.phi,self.mesh.x,self.mesh.y,self.mesh.z,self.ng);
                 self = self.accelerate(); 
              end
              phiold = self.phi;
            end
            
          elseif self.MethodEigs == 1
            self.k = eigs(self.M\self.F);
          end
          
          self.print();
          
        end
        
        function self=checkParams(self)
            %%%%%%%%%%%%%%%%Build Problem Geometry%%%%%%%%%%%%%%%%%%%%
            %Evenly space mesh points
            if (self.dimensions > 3 || self.dimensions < 1)
                error('Invalid number of dimensions');
                
            end

            % if (dimensions < 3)
            %    if (z_len ~= 1)
            %        disp('Invalid number of mesh cells in z-axis.');
            %        disp('Less than 3D. Number of cells in z-axis must...
            %          be set to 1');
            %        break
            %    end
            %    dim(3) = 1;
            %    BC(3) = 1;
            % end
            % 
            % if (dimensions < 2)
            %    if (y_len ~= 1)
            %        disp('Invalid number of mesh cells in y-axis.');
            %        disp('Problem is 1D. Number of cells in y-axis must...
            %          be set to 1');
            %        break
            %    end
            %    dim(2) = 1;
            %    BC(2) = 1;
            % end
            self.mesh = self.mesh.setEqualDist(1,self.dim(1),self.dim(2),self.dim(3));
            %Set boundary conditions
            self.mesh = self.mesh.setBC(self.BC(1), self.BC(2), self.BC(3));%set BC
            %Compute total number of mesh cells
            self.total_mesh = self.mesh.x * self.mesh.y * ...
                self.mesh.z * self.mesh.g;
            
        end
        
        
        function self=accelerate(self)
          
          self = self.calculateCurrents();
          
          %Make coarse mesh from fine mesh
          coarse = Mesh.fineToCoarse(self.mesh,self.phi,self.gridReductionFactor, self.dimensions);
          
          coarse_flux_old = coarse.phi;
          
          [self.coarseM coarse] = BuildLossMatrix(coarse, 1, true, self.mesh, 1.0);
          self.coarseF = BuildProductionMatrix(coarse, 1);

          total_mesh_coarse = coarse.x * coarse.y * coarse.z * coarse.g;
          %Initial Guess:
          format long
          k2 = 1;
          kold2 = .1;
          phiold = ones(total_mesh_coarse,1);
          for iters_coarse = 1:self.max_iters_coarse

              
              if (self.innerprintlevel == 1)
                  text = strcat('CMFD Iteration: ',num2str(iters_coarse),'. Error in k: ',num2str(abs(k2-kold2)));
                  disp(text);
              end
              if (abs(k2-kold2) < self.convergence_coarse)
                  break;
              end
              kold2 = k2;
              b2 = self.coarseF*phiold/kold2;
              phi2 = self.coarseM\b2;
              k2 = kold2 * sum(phi2)/sum(phiold);
              phi2 = phi2 / sum(sum(phi2));
              phiold = phi2;
          end
          
          %CMFD completed. Remap back onto fine mesh
          %self.phi = reshape(self.phi,self.mesh.x,self.mesh.y,self.mesh.z,self.ng);
          %phi2 = reshape(phi2,coarse.x,coarse.y,coarse.z,coarse.g);
          info = [self.dimensions self.mesh.g self.mesh.x self.mesh.y];
          info2 = [self.dimensions coarse.g coarse.x coarse.y];
            for i = 1:self.mesh.x
                for j = 1:self.mesh.y
                    for k = 1:self.mesh.z
                        for g = 1:self.ng
                        
                            %find coarse mesh indices
                            i2 = ceil(i/self.gridReductionFactor(1));
                            j2 = ceil(j/self.gridReductionFactor(2));
                            k2 = ceil(k/self.gridReductionFactor(3));
                            g2 = g;
                            
            self.phi(indexToMat(i,j,k,g,info)) = self.phi(indexToMat(i,j,k,g,info)) *...
                phi2(indexToMat(i2,j2,k2,g2,info2)) / coarse_flux_old(i2,j2,k2,g2);
                %coarse.phi(i2,j2,k2,g2) / coarse_flux_old(i2,j2,k2,g2);           
                        end
                    end
                end
            end    
          self.phi = reshape(self.phi,self.total_mesh,1);
          
        end
        
        function self = calculateCurrents(self) 
            
            INF = [self.dimensions self.mesh.g self.mesh.x self.mesh.y];
            %calculate all surface currents
            for i1 = 1:self.mesh.x
                for j1 = 1:self.mesh.y
                    for k1 = 1:self.mesh.z
                        for g1 = 1:self.ng
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
                                neigh = checkBoundary(i1,j1,k1,self.mesh.x,self.mesh.y,self.mesh.z,xyz,side);
                                if (neigh == 1)
                                %not boundary
                                    %calculate current for each boundary
                                    self.mesh.Jsurf(i1,j1,k1,g1,surf) = -self.mesh.dtilde(i1,j1,k1,g1,surf) * ...
                                        (self.phi(indexToMat(i1,j1,k1,g1,INF)) - self.phi(indexToMat(i2,j2,k2,g1,INF)));                                
                                else
                                %boundary
                                %Only works for vacuum currently
                                    self.mesh.Jsurf(i1,j1,k1,g1,surf) = -self.mesh.dtilde(i1,j1,k1,g1,surf) * ...
                                        self.phi(indexToMat(i1,j1,k1,g1,INF));
                                end
                                if (side == 1)
                                        self.mesh.Jsurf(i1,j1,k1,g1,surf) = -self.mesh.Jsurf(i1,j1,k1,g1,surf);
                                end
                            end
                        end
                    end
                end
            end    
            
        end
        
        function self = print(self)
            
            disp k;
            disp (self.k);
            disp('iters');
            disp(self.iters);

            if (self.MethodGauss && self.dimensions == 1)
                plotflux = reshape(self.phi,self.ng,numel(self.phi)/self.ng);
                plot(self.mesh.X,plotflux(1,:));
                hold on
                if (self.ng == 2) 
                    plot(self.mesh.X,plotflux(2,:),'r');
                end
                xlabel('x (cm)')
                ylabel('Normalized Flux');
                legend('Group 1', 'Group 2');
                title('Reactor Flux Profile');
            end

            %if (self.MethodGauss && self.dimensions == 1 && self.CMFD == 1)
            %    flux2 = reshape(phi2,self.ng,numel(phi2)/self.ng);
            %    figure;
            %    plot(coarse.X,flux2(1,:));
            %    hold on
            %    if (self.ng == 2) 
            %        plot(coarse.X,flux2(2,:),'r');
            %    end
            %    xlabel('x (cm)')
            %    ylabel('Normalized Flux');
            %    legend('Group 1', 'Group 2');
            %    title('Reactor Flux Profile - Coarse Mesh');
            %end
            
        end
        
    end
    
  
    
end
	
