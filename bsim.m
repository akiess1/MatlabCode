% main simulation script
% created: 4/20/17 2:28PM
%------------------------------------------------------------------------

function outdata = bsim(xmin,xmax,Nx,ymin,ymax,Ny,...
        tmin,tmax,dt,veldelay,mutdelay,repdelay,...
        difb,dif1,dif2,sec1,sec2,alp1,alp2,beta,sat1,sat2,lam1,lam2,mprob,mdiff,...
        Nb,numGrps,grpsep,aligngrps,sticky,policing,...
        velmax,rotation,vrad,veltype,...
        saveRate,graphics,saveVid,vidName)

% discretize space:
dx = (xmax - xmin)/Nx;
dy = (ymax - ymin)/Ny;
x = (xmin-dx):dx:(xmax+dx); % including ghost nodes x(1) and x(Nx+3)
y = (ymin-dx):dy:(ymax+dx); % including ghost nodes y(1) and y(Ny+3)

Nt = (tmax - tmin)/dt; % number of timesteps
t = tmin; % initialize time

% save variable:
nsave = ceil(Nt/saveRate);
dataRec = cell(1,nsave);

% generate velocity:
[velX,velY] = genVelocity(veltype,velmax,rotation,vrad,x,y);

% initialize fields and bacteria:
bacteria = constructBacteria(xmin,xmax,ymin,ymax,sec1,sec2,...
    Nb,numGrps,grpsep,aligngrps); %*** add policing to this, as well as sticky
% chemicals
c1 = zeros(Nx+3,Ny+3);
c2 = zeros(Nx+3,Ny+3); % including ghost node)

% graphics variables:
if saveVid == 1
    vid = VideoWriter(vidName,'MPEG-4');
    open(vid);
end

ms = max(sec1,sec2);
clrscl = 0.1/(dx*dy*ms*dt);

% figure for ``video''
if graphics == 1
    scrsz = get(groot,'ScreenSize');
    % scrsz = [1 1 width height]
    % position: [left bottom width height]
    h = figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/2]);
end

%-------------------------------------------------------------------------
% LOOP THROUGH TIME:
alive = 1;
for ti = 1:Nt
    % test delay:
    if t >= veldelay
        delaycntrl = 1;
    else
        delaycntrl = 0;
    end
    if t >= mutdelay
        mdel = 1;
    else 
        mdel = 0;
    end
    % BACTERIA DYNAMICS:
    [cp1, cp2] = bsecrete(xmin,xmax,Nx,ymin,ymax,Ny,bacteria); 
    bacteria = moveBacteria(xmin,xmax,Nx,ymin,ymax,Ny,dt,...
        difb,velX,velY,bacteria,...
        veltype,rotation,vrad,delaycntrl,sticky);
    if ti > repdelay
        %fprintf('time = %d \n',ti);
        [bacteria,alive] = reproduceBacteria(xmin,xmax,Nx,ymin,ymax,Ny,dt,...
            alp1,alp2,beta,sat1,sat2,mdel*mprob,mdiff,bacteria,c1,c2);
    end % adjust this later -- reproduction delay
 %   fprintf('nb = %d \n',length(bacteria)); % as a check***
    if alive == 0
        disp('Everybody died!');
        break;
    end
    
    % CHEMICAL DYNAMICS:
    % boundary conditions:
    c1(1,:) = c1(Nx+2,:); % periodic
    c1(Nx+3,:) = c1(2,:); % periodic
    c1(:,1) = c1(:,2); % no flux
    c1(:,Ny+3) = c1(:,Ny+2); % no flux
    
    c2(1,:) = c2(Nx+2,:); % periodic
    c2(Nx+3,:) = c2(2,:); % periodic
    c2(:,1) = c2(:,2); % no flux
    c2(:,Ny+3) = c2(:,Ny+2); % no flux
    
    % chem 1:
    cppdiff1 = chemDiffusion(xmin,xmax,Nx,ymin,ymax,Ny,dif1,c1);
    cppadv1 = chemAdvect(xmin,xmax,Nx,ymin,ymax,Ny,...
        delaycntrl*velX,delaycntrl*velY,c1);
    c1 = c1 + (cp1 + cppdiff1 + cppadv1 - lam1*c1)*dt;
    % chem 2:
    cppdiff2 = chemDiffusion(xmin,xmax,Nx,ymin,ymax,Ny,dif2,c2);
    cppadv2 = chemAdvect(xmin,xmax,Nx,ymin,ymax,Ny,...
        delaycntrl*velX,delaycntrl*velY,c2);
    c2 = c2 + (cp2 + cppdiff2 + cppadv2 - lam2*c2)*dt;
    
    % update time:
    t = t + dt;
    
    % save data:
    dataRec{ti} = bacteria;
    
    % display graphics:
    if graphics == 1
        displayGraphics(xmin,xmax,Nx,ymin,ymax,Ny,bacteria,alive,c1,c2,...
            clrscl,t,h);
    end
   
end % end time loop
%----------------------------------------------

% close video file:
if saveVid == 1
    close(vid);
end

% out put save variable:
outdata = dataRec;
end % END OF FUNCTION





