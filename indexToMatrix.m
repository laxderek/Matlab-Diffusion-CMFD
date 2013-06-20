function [index] = indexToMatrix(dim,i,j,k,g,ng,nx,ny)

    if (dim == 1)
        index = (g - 1) + ng*(i - 1);
        if (j ~= 1 || k ~= 1)
            index = 0;
        end
    elseif (dim == 2)
        index = (g - 1) + ng*(i - 1) + ng*nx*(j - 1) + 1;
        if (k ~= 1)
            index = 0;
        end
    else
        index = (g - 1) + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1) + 1;
    end
    
    if (g < 1 || g > ng)
        index = 0;
    end
    if (i < 1 || i > nx)
        index = 0;
    end
    if (j < 1 || j > ny)
        index = 0;
    end
    if (k < 1)
        index = 0;
    end
end