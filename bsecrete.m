% secrete onto chemical grid:

function [cp1, cp2] = bsecrete(xmin,xmax,Nx,ymin,ymax,Ny,bacteria)

ctmp1 = zeros(Nx+3,Ny+3);
ctmp2 = zeros(Nx+3,Ny+3); % including ghost node

dx = (xmax - xmin)/Nx;
dy = (ymax - ymin)/Ny;
invArea = (dx*dy); % ***double check this...

Nb = length(bacteria);

for bi = 1:Nb
    xp = bacteria(bi).position(1);
    yp = bacteria(bi).position(2);
    s1 = bacteria(bi).secretion(1);
    s2 = bacteria(bi).secretion(2);
    [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp);
    ctmp1(xi,yi) = ctmp1(xi,yi) + s1*invArea;
    ctmp2(xi,yi) = ctmp2(xi,yi) + s2*invArea;
end

cp1 = ctmp1;
cp2 = ctmp2;
end