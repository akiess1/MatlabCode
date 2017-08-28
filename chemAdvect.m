   
% ***try higher order upwind schemes later (need more ghosts)...
function cpp = chemAdvect(xmin,xmax,Nx,ymin,ymax,Ny,velX,velY,c)
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;

    mmx = (0.5/dx);
    mmy = (0.5/dy);
    
    % including ghosts there are N+3 nodes along each dimension
    ctmp = zeros(Nx+3,Ny+3);
    order = 1; % first order upwind
    
    if order == 0
    for xi = 2:(Nx+2)
        for yi = 2:(Ny+2)
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi+1,yi)-c(xi-1,yi))*mmx + ...
                    -velY(xi,yi)*(c(xi+1,yi)-c(xi,yi-1))*mmy;
        end
    end
    elseif order == 2
 % possible upwind schemes:   
    
for xi = 3:(Nx+1) % need to add more ghost nodes or something
    for yi = 3:(Ny+1)
        if velX(xi,yi) >= 0 && velY(xi,yi) >= 0
        %if velX(xi,yi) <= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(3*c(xi,yi)-4*c(xi-1,yi)+c(xi-2,yi))*mmx + ...
                -velY(xi,yi)*(3*c(xi,yi)-4*c(xi,yi-1)+c(xi,yi-2))*mmy;
        elseif velX(xi,yi) >= 0 && velY(xi,yi) <= 0
        %elseif velX(xi,yi) <= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(3*c(xi,yi)-4*c(xi-1,yi)+c(xi-2,yi))*mmx + ...
               -velY(xi,yi)*(-c(xi,yi+2)+4*c(xi,yi+1)-3*c(xi,yi))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) >= 0
        %elseif velX(xi,yi) >= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(-c(xi+2,yi)+4*c(xi+1,yi)-3*c(xi,yi))*mmx + ...
                -velY(xi,yi)*(3*c(xi,yi)-4*c(xi,yi-1)+c(xi,yi-2))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) <= 0
       % elseif velX(xi,yi) >= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(-c(xi+2,yi)+4*c(xi+1,yi)-3*c(xi,yi))*mmx + ...
                -velY(xi,yi)*(-c(xi,yi+2)+4*c(xi,yi+1)-3*c(xi,yi))*mmy;
        end
    end
end    
    elseif order == 1
for xi = 2:(Nx+2)
    for yi = 2:(Ny+2)
        if velX(xi,yi) >= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi,yi)-c(xi-1,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi)-c(xi,yi-1))*mmy;
        elseif velX(xi,yi) >= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi,yi)-c(xi-1,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi+1)-c(xi,yi))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) >= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi+1,yi)-c(xi,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi)-c(xi,yi-1))*mmy;
        elseif velX(xi,yi) <= 0 && velY(xi,yi) <= 0
            ctmp(xi,yi) = -velX(xi,yi)*(c(xi+1,yi)-c(xi,yi))*mmx + ...
                -velY(xi,yi)*(c(xi,yi+1)-c(xi,yi))*mmy;
        end
    end
end    
    end
    
cpp = ctmp;
end

