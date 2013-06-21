function [ F ] = BuildProductionMatrix(mesh,dim)

    ng = mesh.g;
    nx = mesh.x;
    ny = mesh.y;
    nz = mesh.z;
    total_mesh = ng * nx * ny * nz;

    %approximated preallocation for speed
    %same for F matrix
    rowF = zeros(total_mesh,1);
    colF = zeros(total_mesh,1);
    valF = zeros(total_mesh,1);

    counterF = 1;
    for irow = 1:total_mesh
        %calc cell specific data
        g = mod(irow-1,ng) + 1;
        i = floor(mod(irow-1,(nx*ng))/ng) + 1;
        if (dim > 1) 
            j = floor(mod(irow-1,(ng*nx*ny)/(ng*nx)) + 1);
        else
            j = 1;
        end
        if (dim > 2) 
            k = floor(mod(irow-1,(ng*nx*ny*nz)/(ng*nx*ny)) + 1);
        else
            k = 1;
        end
            %volume = prod(mesh.dxyz(i,j,k));
        for h = 1:ng
            
            neigh = (h + ng*(i-1) + ng*nx*(j-1) + ng*nx*ny*(k-1));%neighbor matrix index
            nuSigf = mesh.nusigF(i,j,k,g,h);%*volume;
            colF(counterF) = irow;
            rowF(counterF) = neigh;
            valF(counterF) = nuSigf;
            counterF = counterF + 1;     
        end
    end
    rowF(counterF) = total_mesh;
    colF(counterF) = total_mesh;
    valF(counterF) = 0;
    F = sparse(rowF,colF,valF);
end
        
		
        
