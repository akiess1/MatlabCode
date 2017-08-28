 

function fit = findFitness(xi,yi,c1,c2,alp1,alp2,beta,sat1,sat2,s1)

fit = alp1*c1(xi,yi)/(c1(xi,yi)+sat1) -...
        alp2*c2(xi,yi)/(c2(xi,yi)+sat2) - beta*s1;
end

