
% ***rewrite a laplacian function later...
function cpp = chemDiffusion(xmin,xmax,Nx,ymin,ymax,Ny,diff,c)

    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;

    mmx = (1/dx)^2;
    mmy = (1/dy)^2;
    
    % including ghosts there are N+3 nodes along each dimension
    ctmp = zeros(Nx+3,Ny+3);
    
   for xi=2:(Nx+2)
        for yi=2:(Ny+2)
            ctmp(xi,yi) = (c(xi+1,yi) - 2*c(xi,yi) + c(xi-1,yi))*mmx + ...
                (c(xi,yi+1) - 2*c(xi,yi) + c(xi,yi-1))*mmy; 
        end % end if near channels
    end % end loops through r and z
    
  
cpp = diff*ctmp;
end