% print program summary and check if parameters are valid:

function pcheck = printprogsum(xmin,xmax,Nx,ymin,ymax,Ny,...
    tmin,tmax,dt,veldelay,mutdelay,repdelay,...
    difb,dif1,dif2,sec1,sec2,alp1,alp2,beta,sat1,sat2,lam1,lam2,mprob,mdiff,...
    Nb,numGrps,grpsep,aligngrps,sticky,policing,...
    velmax,rotation,vrad,veltype,...
    saveRate,graphics,saveVid,filename,datechar)

fprintf('Program run date: %s \n',datechar);
fprintf('Domain size: [%.2f, %.2f] X [%.2f, %.2f] \n',xmin,xmax,ymin,ymax);
fprintf('Discritization: Nx = %d | Ny = %d \n',Nx,Ny);
fprintf('Run Time: [%.2f, %.2f] | dt = %.6f \n',tmin,tmax,dt);
fprintf('\n*** system parameters *** \n');
fprintf(' db = %.3f | d1 = %.3f | d2 = %.3f \n',difb,dif1,dif2);
fprintf(' s1 = %.3f | s2 = %.3f \n',sec1,sec2);
fprintf(' l1 = %.3f | l2 = %.3f \n',lam1,lam2);
fprintf(' k1 = %.3f | k2 = %.3f \n',sat1,sat2);
fprintf(' a1 = %.3f | a2 = %.3f | beta = %.3f \n',alp1,alp2,beta);
fprintf(' mprob = %.3f | mdiff = %.3f \n',mprob,mdiff);
% veltype options:
% 1 = couette
% 2 = poisuille
% 3 = rankine
% 4 = constant (useful for sticky bacteria)
if veltype == 1
    fprintf(' Veltype: Couette: vmax = %.3f \n',velmax);
elseif veltype == 2
    fprintf(' Veltype: Poiseuille: vmax = %.3f \n',velmax);
elseif veltype == 3
    fprintf(' Veltype: Rankine: rot = %.3f | rad = %.3f \n',rotation,vrad);
elseif veltype == 4
    fprintf(' Veltype: Constant: vmax = %.3f \n',velmax);
else
    fprintf(' no flow velocity ');
end
fprintf(' Initial bnum = %d | initial numgrps = %d \n',Nb,numGrps);
if aligngrps == 1
    fprintf(' groups aligned with separation distance = %.3f \n',grpsep);
end
if sticky == 1
    fprintf(' Using sticky bacteria \n');
end
fprintf(' Velocity delay time = %.5f \n',veldelay);
fprintf(' Mutation delay time = %.5f \n',mutdelay);
fprintf(' Reproduction delay time = %.5f \n',repdelay*dt);
fprintf(' saverate = %d | graphics = %d | video = %d \n',saveRate,graphics,saveVid);
fprintf('Save file name: "%s" \n************************************\n\n',filename);

tmpcheck = 1;
if xmin > xmax || ymin > ymax || tmin > tmax
    tmpcheck = -1;
    disp('bad domain');
end
    
pcheck = tmpcheck;
end