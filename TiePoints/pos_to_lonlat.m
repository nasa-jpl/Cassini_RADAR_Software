function [lon,lat]=pos_to_lonlat(xb);
% function [lon,lat]=pos_to_lonlat(xb);
% positions should be 3 by 1 column vectors  
% or 3 by N matrices of N different positions
% return 1 by N lon and 1 by N lat
uxb=xb./(ones(3,1)*sqrt(sum(xb.*xb)));
theta = acos(uxb(3,:));
lon = atan2(uxb(2,:),uxb(1,:))*180/pi;
lat= pi/2 -theta;
lat= lat*180/pi;
lon=360-lon;
i=find(lon>180); 
lon(i)=lon(i)-360;
i=find(lon<-180); 
lon(i)=lon(i)+360;

