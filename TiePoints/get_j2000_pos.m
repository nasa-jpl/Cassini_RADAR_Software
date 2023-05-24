function [xi]=get_j2000_pos(rx,ry,rz,M,w,dop,range,scpos,scvel,ground,lambda)
%function [xi]=get_j2000_pos(rx,ry,rz,M,w,dop,range,scpos,scvel,ground,lambda)
% routine to determine inertial coordinates of a tiepoint from raw measurement data


% convert input radtio rate w from deg/day to rad/s
dtr=pi/180;
w=(w*dtr)/(24*3600);  % converts degrees_per_day to rads_per_sec

% rotate spacecraft position into Titan body fixed coordinates
% using rotation matrix M
scpos0 = scpos;
scpos=M*scpos;


% Convert s/c velocity into Titan body fixed coordinates
% accounting for the rotation rate of Titan
% As it happens this is insensitive to the small spin rate variations
% within the realm of possible Titan spin models
scvel=M*scvel+[0,+w,0;-w,0,0;0,0,0]*scpos;
speed=norm(scvel);

% compute dot product of velocity and look vector from Doppler
vdotr=dop*lambda/2.0;

% compute two points in body fixed coordinates where the range sphere, Titan
% surface sphere (adjusted for SARTopo) and Doppler cone intersect
points=ellsphcone_intersect(scpos,scvel,vdotr,range,rx,ry,rz);


sz=size(points);
N=sz(2);

mindist=10000;

for i=1:N

% rotate each position back to inertial coordinates
x=(M^-1)*points(:,i);


%%%%% choose point  closest ground location unit vector --- ground
%%%%% the ground vector used for the purpose doesn't need to be accurate
if(isreal(x))
   dist=norm(x/norm(x)-ground);
else
   dist = 100000;
end
if(dist<mindist)
     xi=x;
     mindist=dist;
end
end
