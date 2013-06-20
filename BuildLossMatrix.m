function [M mesh] = BuildLossMatrix(mesh, dim, CMFD, fine, damping)
	

    ng = mesh.g;
    nx = mesh.x;
    ny = mesh.y;
    nz = mesh.z;
    total_mesh = mesh.x * mesh.y * mesh.z * mesh.g;

    
    %matrices to turn into sparse matrix for M
    %approximated preallocation for speed
    rowM = zeros(total_mesh*2,1);
    colM = zeros(total_mesh*2,1);
    valM = zeros(total_mesh*2,1);
    counterM = 1;
    %loop over rows
    for irow=1:total_mesh  %indexed by
    %for irow=1:1
        %for each row, there are a number of contributions to the cell
        %these will be compiled in 3 lists and used to build a sparse
        %matrix at the end


        %for neighbors
        %addrow = 
        %row = [row addrow];
        %col = [col addcol];
        %val = [val addval];

        g = mod(irow - 1, ng) + 1; %current group number
        i = floor(mod(irow - 1, ng*nx)/ng) + 1; % current x loc

        if (dim > 1)
            j = floor(mod(irow - 1, ng*nx*ny)/(ng*nx))+ 1;
        else 
            j = 1;
        end
        if (dim > 2)
            k = floor(mod(irow - 1, ng*nx*ny*nz)/(ng*nx*ny)) + 1;
        else 
            k = 1;
        end
        %pull material properties from Mesh object based on physical location
        %volume = prod(mesh.dxyz(i,j,k));
        %index = (int)FastMath.ceil((g+ng*i + ng*nx*j + ng*nx*ny*k)/2);
        totalxs = mesh.sigT(i,j,k,g);
        scatxs = mesh.sigS(i,j,k,g,g);
        diffusion = mesh.D(i,j,k,g);

        dhat = 0; %%%%%%%ENABLE FOR CMFD

        for surf = 1:6          % 1 = -x; 2 = +x
                                % 3 = -y; 4 = +y
                                % 5 = -z; 6 = +z
                                %LEAKAGE

            %%determine x, y, or z surface
            xyz = round((surf / 2)); %1=x;2=y;3=z
            %%determine left or right surface
            side = -2 * mod(surf,2) + 1; % -1=left, 1=right
            %%determine neighbor index
            %if doesn't exist (boundary), returns -1
            if (xyz == 1)
                i2 = i + side;
                j2 = j;
                k2 = k;
            elseif (xyz == 2)
                i2 = i;
                j2 = j + side;
                k2 = k;
            else %(xyz == 3) 
                i2 = i;
                j2 = j;
                k2 = k + side;		        	
            end
            neigh = getNeighbor(irow,ng,nx,ny,total_mesh, xyz, side);
            if (neigh == -1) %boundary
                %apply boundary condition            
                mesh.dtilde(i,j,k,g,surf)= 2 * diffusion * (1 - mesh.albedo(xyz)) / ...
                    (4*diffusion*(1 + mesh.albedo(xyz)) + mesh.dxyz(xyz) * (1 - mesh.albedo(xyz)));
                if CMFD == true
                   %calculate dhat
                   if (i == 1) %left boundary
                   i3 = i*mesh.gridreduce(1) - 1;
                   else %right boundary
                   i3 = i*mesh.gridreduce(1);                       
                   end
                   j3 = j*mesh.gridreduce(1);                   
                   k3 = k*mesh.gridreduce(1);
                   if (dim < 3)
                       k3 = k;
                   end
                   if (dim < 2)
                       j3 = j;
                   end
                   dhatnew = -(fine.Jsurf(i3,j3,k3,g,surf) - side *  mesh.dtilde(i,j,k,g,surf)* ...
                       mesh.phi(i,j,k,g)) / mesh.phi(i,j,k,g);
                   dhatold = mesh.dhat(i,j,k,g,surf);
                   
                   dhat = dhatnew * damping + (1 - damping) * dhatold;
                   mesh.dhat(i,j,k,g,surf) = dhat;
                end
            else %not boundary
                %apply current
                rowM(counterM) = irow; 
                colM(counterM) = neigh;
                neigh_diff = mesh.D(i2,j2,k2,g);
                mesh.dtilde(i,j,k,g,surf) = 2 * diffusion * neigh_diff / ...
                    (neigh_diff*mesh.dxyz(xyz) + diffusion*mesh.dxyz(i2,j2,k2,xyz));
                if CMFD == true
                   %calculate dhat
                        i3 = i*mesh.gridreduce(1);
                        j3 = j*mesh.gridreduce(2);
                        k3 = k*mesh.gridreduce(3);
                    if (xyz == 1 && side == -1)
                        i3 = (i+side)*mesh.gridreduce(1) + 1;
                    elseif (xyz == 2 && side == -1)
                        j3 = (j+side)*mesh.gridreduce(2) + 1;
                    elseif (xyz == 3 && side == -1) 
                        k3 = (k+side)*mesh.gridreduce(3) + 1;
                    end
                   if (dim < 3)
                       k3 = k;
                   end
                   if (dim < 2)
                       j3 = j;
                   end
                   current = fine.Jsurf(i3,j3,k3,g,surf);
                   if (surf == 1 || surf == 2)                  
                   dhatnew = -(current + side * mesh.dtilde(i,j,k,g,surf) * ...
                       (-mesh.phi(i,j,k,g) + mesh.phi(i2,j2,k2,g))...
                       )/ (mesh.phi(i,j,k,g) + ...
                       mesh.phi(i2,j2,k2,g));
                   dhatold = mesh.dhat(i,j,k,g,surf);
                   
                   dhat = dhatnew * damping + (1 - damping) * dhatold;
                   mesh.dhat(i,j,k,g,surf) = dhat;
                   end
                end
                
                
                valM(counterM) = (-mesh.dtilde(i,j,k,g,surf) + side * dhat) / mesh.dxyz(i,j,k,xyz);
                counterM = counterM + 1;
            end 
            mesh.Leakage(i,j,k,g,surf) = mesh.dtilde(i,j,k,g,surf) - side * dhat;
        end

        jtot = (mesh.Leakage(i,j,k,g,2) + mesh.Leakage(i,j,k,g,1)) / mesh.dxyz(i,j,k,1) + ...
                  (mesh.Leakage(i,j,k,g,4) + mesh.Leakage(i,j,k,g,3)) / mesh.dxyz(i,j,k,2) + ...
                  (mesh.Leakage(i,j,k,g,6) + mesh.Leakage(i,j,k,g,5)) / mesh.dxyz(i,j,k,3);
        %sum(mesh.Leakage(i,j,k,g,1:6));

        colM(counterM) = irow;
        rowM(counterM) = irow;
        valM(counterM) = jtot + (totalxs - scatxs);
        counterM = counterM + 1;
        
        for h = 1:ng        

            neigh = (h + ng*(i-1) + ng*nx*(j-1) + ng*nx*ny*(k-1));%neighbor matrix index
            %if in-scattering (group = g), do nothing
            if (g == h)
                continue;
            else
            %do group -> group scattering
            rowM(counterM) = irow;
            colM(counterM) = neigh;
            valM(counterM) = -mesh.sigS(i,j,k,h,g);
            counterM = counterM + 1;
            end        

        end
   
    
    end
    
        M = sparse(rowM,colM,valM);
end
