function [index] = indexToMat(i,j,k,g,info)
    
    if (info(1) == 1)
        index = (g - 1) + info(2)*(i - 1) + 1;
        if (j ~= 1 || k ~= 1)
            index = 0;
        end
    elseif (info(1) == 2)
        index = (g - 1) + info(2)*(i - 1) + info(2)*info(3)*(j - 1) + 1;
        if (k ~= 1)
            index = 0;
        end
    else
        index = (g - 1) + info(2)*(i - 1) + info(2)*info(3)*(j - 1) + info(2)*info(3)*info(4)*(k - 1) + 1;
    end
    
    if (g < 1 || g > info(2))
        index = 0;
    end
    if (i < 1 || i > info(3))
        index = 0;
    end
    if (j < 1 || j > info(4))
        index = 0;
    end
    if (k < 1)
        index = 0;
    end
end