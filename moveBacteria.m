   

function bact = moveBacteria(xmin,xmax,Nx,ymin,ymax,Ny,dt,...
        difb,velX,velY,bacteria,veltype,rotation,vrad,...
        delay,sticky)
    
    Nb = length(bacteria);
    dx = (xmax - xmin)/Nx;
    dy = (ymax - ymin)/Ny;

    for bi = 1:Nb
        xp = bacteria(bi).position(1);
        yp = bacteria(bi).position(2);
        x = [xp yp];
        if veltype == 3
            % get 'exact' velocity for rankine vortex
            [bvx, bvy] = brankine(rotation,vrad,xp,yp);
        else
            %*** maybe get exact velcity in general?
            [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp);
            bvx = velX(xi,yi);
            bvy = velY(xi,yi); 
        end
        if sticky ==  1
            vel = [0 0];
        else
            vel = delay*[bvx bvy];  
        end
        
        probwalk = randi(4);
        stepwalk = [0 0];
        if probwalk == 1
            stepwalk = [1 0];
        elseif probwalk == 2
            stepwalk = [-1 0];
        elseif probwalk == 3
            stepwalk = [0 1];
        elseif probwalk == 4
            stepwalk = [0 -1];
        end
        dpos = (difb*stepwalk*(dt/(dx*dy))) + (vel*dt);
        
        % update position:
        newPos = x + dpos;
        newX = newPos(1);
        newY = newPos(2);
    
        % check boundary conditions:
        if newX < xmin
            % periodic bc:
            xdif = newX - xmin; % this is negative
            newX = xmax + xdif; % this is less than xmax
        elseif newX > xmax
            xdif = newX - xmax; % this is positive
            newX = xmin + xdif; % this is larger than xmin
        end % possible problems when movement is larger than domain 
        %(but that shouldn't occur)
        if newY < ymin
            % reflective bc:
            ydif = ymin - newY; % this is positive
            newY = ymin + ydif; % this is larger than ymin
        elseif newY > ymax
            ydif = newY - ymax; % this is positive
            newY = ymax - ydif; % this is less than ymax
        end
        % fix position:
        bacteria(bi).position(1) = newX;
        bacteria(bi).position(2) = newY;
    end

bact = bacteria;
end