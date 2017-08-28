 

 %  ***later add mutation of diffusion constant 
  function  [bact,state] = reproduceBacteria(xmin,xmax,Nx,ymin,ymax,Ny,dt,...
        alp1,alp2,beta,sat1,sat2,mprob,mdiff,bacteria,c1,c2)
    
    Nb = length(bacteria);
    numRep = 0;
    toRep = zeros(2,Nb); % Nb array of vectors: [index, number to reproduce]
    alive = 1;
    
    k = 1;
    for bi = 1:Nb
        xp = bacteria(bi).position(1);
        yp = bacteria(bi).position(2);
        sec1 = bacteria(bi).secretion(1);
        [xi,yi] = findGridIndex(xmin,xmax,Nx,ymin,ymax,Ny,xp,yp);
        fit = findFitness(xi,yi,c1,c2,alp1,alp2,beta,sat1,sat2,sec1)*dt;

        fitInt = floor(fit);
        fitDec = fit - fitInt;

        % reproduce with probability prepr 
        prepr = rand;
        fitProb = 0;
        if prepr < fitDec
            fitProb = 1;
        end
        if fitInt >= 0
            nrep = fitInt + fitProb + 1;
        else
            nrep = 0;
        end % if fitness is positive

        % fit = number to reproduce = fitInt + fitProb

        if fit > 0
          toRep(1,k) = bi;
          toRep(2,k) = nrep;
          numRep = numRep + nrep;
          k = k + 1;
        end
    end % for each bacteria -- record number of offspring to produce

    if numRep == 0
        alive = 0;
        state = alive;
        bact = 0;
        return;
    end %if bacteria go extinct
    
    % === contruct array of new bacteria: ===
    temp(numRep) = struct('position', zeros(1,2), ...
                        'secretion', zeros(1,2), ...
                        'mutations',0);

    ri = 1; % index for toRep(1,:) array
    ti = 1; % index for temp(:) array

    Nreproducing = length(toRep(1,:));

    ind = toRep(1,ri);
    numOffspring = toRep(2,ri);
    while ind ~= 0
        for i = 1:numOffspring
            temp(ti).position = bacteria(ind).position;
            temp(ti).secretion = bacteria(ind).secretion;
            temp(ti).mutations = bacteria(ind).mutations;
            %mutate offspring with probability PrMt:
             prob = rand;
             if prob < mprob
                sctp = temp(ti).secretion(1);
                sgn = sign(sctp);
                if sgn == 0
                 sgn = 1;
                end % temporary solution to evolution of zero (always + for sinpop)
                dsc = mdiff*(2*rand - 1);
                sctp = sctp + dsc;
                if sign(sctp) ~= sgn
                 sctp = -sctp;
                end
                temp(ti).secretion(1) = sctp;
                temp(ti).mutations = temp(ti).mutations + 1;
             end % mutate with probability PrMt
            % increment temp index:
             ti = ti + 1;
        end % for number of reproductions of bacteria `ind'
        % increment producer index:
        ri = ri + 1;
        if ri > Nreproducing 
            break;
        end
        ind = toRep(1,ri);
        numOffspring = toRep(2,ri);
    end % while index refers to bacteria to reproduce

    % is this necessary?:
     NewNb = length(temp);
     % shuffle array for purposes of secretion randomness.....
     shuffleInd = randperm(NewNb,NewNb);
     temp = temp(shuffleInd);

  state = alive;  
  bact = temp;
  end % END OF FUNCTION
  
  