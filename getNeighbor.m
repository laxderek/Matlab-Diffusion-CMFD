function [index] = getNeighbor(irow,ng,nx,ny,total_mesh,xyz,side)

    if xyz == 1
        index = irow + side * ng;
    elseif xyz == 2
        index = irow + side * ng * nx; 
    else
        index = irow + side * ng * nx * ny;
    end

    if index < 1 || index > total_mesh
        index = -1;
    end
    if irow
end