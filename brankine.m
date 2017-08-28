


function [bvx, bvy] = brankine(rotation,radius,xp,yp)

PI = 3.141592653589793238462643383279502884197169399;

dist = sqrt(xp^2 + yp^2);
theta = atan2(xp,yp);

if dist <= radius
	factor = (rotation*dist)/(2*PI*(radius^2));
else
	factor = rotation/(2*PI*dist);
end

bvy = -sin(theta)*factor;
bvx = cos(theta)*factor;
end


