function [g,i,j,k] = matrix_to_indices(irow,ng,nx,ny)

    % compute individual indices
    k = ceil(irow/(ng*nx*ny));
    j = ceil(irow/(ng*nx) - (k-1)*ny);
    i = ceil(irow/ng - (k-1)*ny*nx - (j-1)*nx);
    g = irow - (k-1)*ny*nx*ng - (j-1)*nx*ng - (i-1)*ng;

end