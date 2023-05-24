function [points]=ellsphcone_intersect(scpos,scvel,vdotr,range,rx,ry,rz);
% routine to found two intersection points of 
% 1) sphere with radius rx centered at the origin
% 2) sphere with radius range, centered at scpos
% 3) cone with apex at scpos, center axis along scvel, and vdotr is the
% cosine of half the vertex anglemake_tiep
% (This routine has only been tested for rx=ry=rz)

rs = (scpos(1))^2 + (scpos(2))^2 + (scpos(3))^2;
rs_2 = 0.5*(rx*rx - range*range - rs);
r2 = range * range;
v_r = vdotr * range;

a1 = scpos(1);
a2 = scpos(2);
a3 = scpos(3);
b1 = scvel(1);
b2 = scvel(2);
b3 = scvel(3);
u = rs_2;
v = v_r;

t = a2*b3 - b2*a3;
f1 = (b3*u - a3*v)/t;
f2 = (a3*b1 - a1*b3)/t;
g1 = (a2*v - b2*u)/t;
g2 = (a1*b2 - a2*b1)/t;

d1 = 1 + f2*f2 + g2*g2;
d2 = 2.0*(f1*f2 + g1*g2);
d3 = f1*f1 + g1*g1 - r2;
d = sqrt(d2*d2 - 4.0*d1*d3);
x1 = (-d2 + d)/2.0/d1;
x2 = (-d2 - d)/2.0/d1;

y1 = f1 + f2*x1;
y2 = f1 + f2*x2;
z1 = g1 + g2*x1;
z2 = g1 + g2*x2;

N = 2;

points(1,1) = x1 + scpos(1);
points(2,1) = y1 + scpos(2);
points(3,1) = z1 + scpos(3);

points(1,2) = x2 + scpos(1);
points(2,2) = y2 + scpos(2);
points(3,2) = z2 + scpos(3);

