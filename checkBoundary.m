function [boolVal] = checkBoundary(x,y,z,nx,ny,nz,xyz,side)

    boolVal = 1;
    if (xyz == 1)
        x2 = x + side;
        if (x2 < 1 || x2 > nx)
            boolVal = 0;
        end
    elseif (xyz == 2)
        y2 = y + side;
        if (y2 < 1 || y2 > ny)
            boolVal = 0;
        end
    else %(xyz == 3)
        z2 = z + side;	
        if (z2 < 1 || z2 > nz)
            boolVal = 0;	  
        end
    end



end