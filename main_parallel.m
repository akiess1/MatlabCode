% --------------------------------------------------------------
% MAIN SCRIPT
% AUTHOR: GURDIP S. UPPAL
% DATE: 4/20/17
% LAST MODIFIED: 4/28/17 
% --------------------------------------------------------------
% === CLEAR WORKSPACE ===
close all; clear; clc;
% TIME SIMULATION:
tic; % start stopwatch
%----------------------------------------------------------------

% RUN IN PARALLEL:
parpool(12);
nruns = 450;

parfor rni = 1:nruns

velocities = [0 10 20 30 40 50 60 70 80];

vi = mod(rni-1,9)+1;


% === SYSTEM PARAMETERS ===
a1 = 75;        % resource benefit
a2 = 80;        % waste harm
b1 = 0.2;       % resource secretion cost
k1 = 0.01;      % resource saturation constant
k2 = 0.1;       % waste saturation constant
s1 = 100;       % initial resource secretion (subject to mutations)
s2 = 100;       % waste secretion rate
l1 = 50;        % resource decay rate
l2 = 15;        % waste decay rate
d1 = 5;         % resource diffusion contant
d2 = 15;        % waste diffusion constant
db = 5.0;       % bacteria diffusion (we rescale space to take this out)
% relevant rescaled parameters are given below...

% === SPACE-TIME PARAMETERS ===
% domain:
xmin = -25;
xmax = 25;
ymin = -25;
ymax = 25;
% discretization:
Nx = round((xmax-xmin)*5); %60;   % x nodes - 1
Ny = round((ymax-ymin)*5); %60;   % y nodes - 1
%time:
tmin = 0;
tmax = 0.1; %4.0;
dt = 0.0001; % for stability

% velocity parameters:
rotation = 10000;
vrad = 10;
vmax = velocities(vi);
veltype = 4; % see options below:
% veltype options:
% 1 = couette
% 2 = poisuille
% 3 = rankine
% 4 = constant (useful for sticky bacteria)
% else: no flow - (equivalent to setting vmax = 0)

% mutation parameters:
mutProb = 0.0000000;       % mutation probability; between 0 and 1
mutDiff = 0.1*s1;          % how much secretion rate changes on evolution

% === BACTERIA INITIALIZATION ===
% types: (1 = true)
sticky = 1;        % sticky bacteria do not advect with flow
policing = 0;%*** still working on this - actually policing should be...
% part of the bacteria structure and perhaps it and stickiness can be...
% mutatible qualities
% policing bacteria adjust secretion rate to eliminate cheaters
% number:
numGrps = 1;           % number of initial groups. enter 0 for no groups (random distribution)
Nb = 50*numGrps;        % initial number of bacteria
% align:
aligngrps = 1;          % 1 for aliging groups in row (for sticky sims) else random
grpsep = 2;             % for two aligned groups (for sticky sims)
% dynamics delays:   (to let groups initialize)
veldelay = 0.005;       % delay time for flow 
mutdelay = 0.005;       % delay for mutation 
repdelay = 3; % in time steps (generations)

% === SAVE OPTIONS ===
saveRate = 20;           % collect ccdata modulo every ''saveRate'' time step
graphics = 0;            % 0 = don't display graphics (for crc runs); 1 = display on
saveVid = 0;             % enter 1 to save video after simulation ends
filename = sprintf('sticky_ng%d_v%d_sd%.2f',numGrps,vmax,grpsep);        % name of file to save video
datechar = datestr(now,'mm-dd-yy_HH-MM-SS');
vidName = sprintf('%s_%s_vid.avi',filename,datechar);
dataFile = sprintf('%s_%s_data.mat',filename,datechar);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% === RESCALE PARAMETERS ===
%{
% need to be careful with space and time rescaling...
dif1 = d1/db;
dif2 = d2/db;
sec1 = s1/(k1*l1);
sec2 = s2/(k2*l2);
alp1 = a1/l1;
alp2 = a2/l2;
beta = b1*k1;
sigma = l2/l1;
% mutation rescaling:
mprob = mutProb*((k1*l1)^2);
mdiff = sec1;
% velocity rescaling;
velmax = vmax*sqrt(1/(l1*db)); % or do we take vmax to be the value after rescaling??
%}
% ===== RUN PROGRAM =====

    data = bsim(xmin,xmax,Nx,ymin,ymax,Ny,...
        tmin,tmax,dt,veldelay,mutdelay,repdelay,...
        db,d1,d2,s1,s2,a1,a2,b1,k1,k2,l1,l2,mutProb,mutDiff,...
        Nb,numGrps,grpsep,aligngrps,sticky,policing,...
        vmax,rotation,vrad,veltype,...
        saveRate,graphics,saveVid,vidName);
    parsave(dataFile,data);

end % end of parfor loop


rntime = toc;
fprintf('\nProgram run time: %f \n',rntime);

delete(gcp('nocreate'));
% *END SCRIPT*

