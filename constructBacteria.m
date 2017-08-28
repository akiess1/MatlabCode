% construct bacteria:

function bact = constructBacteria(xmin,xmax,ymin,ymax,sec1,sec2,...
    Nb,numGrps,grpsep,aligngrps)

    % ==== Bacteria Construction ====
    bacteria(Nb) = struct('position', zeros(1,2), ...
                        'secretion', zeros(1,2), ...
                        'mutations',0);

    if numGrps > 0
        if aligngrps == 1
            grpPosY = 0.5*(ymax-ymin)*ones(numGrps,1) + ymin;
            if numGrps == 2
                grpPosX = [0.5*(xmax-xmin-grpsep)+xmin; ...
                    0.5*(xmax-xmin+grpsep)+xmin];
            elseif numGrps == 1
                grpPosX = 0.5*(xmax-xmin)+xmin; % center singleton group
            else
                tmpspacing = linspace(xmin,xmax,numGrps+1);
                grpPosX = transpose(tmpspacing(1:numGrps));
            end
        else
            grpPosX = (xmax-xmin)*rand(numGrps,1) + xmin;
            grpPosY = (ymax-ymin)*rand(numGrps,1) + ymin;
        end % if aligning groups

        grpPos = [grpPosX grpPosY];
    end % create group positions if needed

    for bi = 1:Nb
        if numGrps > 0
            j = mod(bi,numGrps) + 1;
            p1 = grpPos(j,1);
            p2 = grpPos(j,2);
        else
            p1 = (xmax-xmin)*rand + xmin;
            p2 = (ymax-ymin)*rand + ymin;
        end
        pos = [p1 p2];
        bacteria(bi).position = pos;
        bacteria(bi).secretion = [sec1, sec2];
        bacteria(bi).mutations = 0;
    end % for each bacteria

bact = bacteria;
end %END FUNCTION

