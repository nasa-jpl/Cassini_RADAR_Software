function [resid,tiepts,debug]=get_residual_for_case(model,useonlybadties,bruce_outliers)
%%% function [resid,tiepts,debug]=get_residual_for_case(model,useonlybadties,bruce_outliers)
%%% Computes a structure that contains tiepoint lats and lon and residual distance between tiepoint observations
%%% using a spin model structure
%%% Inputs
%%% model: spin model structure
%%% useonlybadties: switch to use good quality tiepoints if set to 0 and bad quality tiepoints if set to 1
%%% bruce_outliers: indices of tiepoints tht failes Bruce Bills' quality control algorithm
%%% Outputs
%%% resid: structure containing residuals distance between tiepoint observations
%%% tiepts: tiepoint information structure
%%% debug: debugging info structure


  
%%% load Bryan Stiles' tiepoint quality metric for each tiepoint
tqual=load('-ASCII','tqual_2016.txt');
tqual_thresh=0.5;

% constants for converting to/from degrees/radians
dtr= pi/180;
rtd= 1/dtr;

% time unit conversions
second_to_day=1/(24*3600);
second_to_cent=1/(100*365.25*24*3600);



%%% read the raw tiepoint data 
[dop1,range1,scpos1,scvel1,dc1,t1,lambda1,dop2,range2,scpos2,scvel2,dc2,t2,lambda2,h2,names]=read_dopran_old('jplusgscornell_2016_newsartopo_randop.tab');
 


%% Compute indices of tiepoints that pass the Stiles and Bills quality controls 

gg=find(tqual>0 & tqual < tqual_thresh);

if(~isempty(bruce_outliers))
  tqual(gg(bruce_outliers))=1+tqual_thresh;
  gg=find(tqual>0 & tqual < tqual_thresh);
end



%% special case when useonlybadties is set to 1
%% resets tiepoint mask to only th bad tiepoints
if(exist('useonlybadties','var'))
  if(useonlybadties)
    gg=find(tqual<=0 | tqual>= tqual_thresh)
  end  
end  

%%% Subset to obtain all the desired good (or bad) tiepoint data
  %%%% Tie point values are dop=doppler, range=range to target, scpos=spacecraft position, scvel = spacecraft velocity, 
  %%%% dc= approximate centroid location of beam,t=time,lambda = wavelength of tranmitted signal
  %%%% h = sartopo height.
  %%%%  1 means first observation 2 means second.
dop1=dop1(gg);
range1=range1(gg);
scpos1=scpos1(:,gg);
scvel1=scvel1(:,gg);
dc1=dc1(:,gg);
t1=t1(gg);
lambda1=lambda1(gg);

dop2=dop2(gg);
range2=range2(gg);
scpos2=scpos2(:,gg);
scvel2=scvel2(:,gg);
dc2=dc2(:,gg);
t2=t2(gg);
lambda2=lambda2(gg);
h2=h2(gg);



%% h1 is the height of the SARTopo from the first flyby
%% h2 is the height from the second. h2 is always used becasue the initial flybies were more errorprone for sartopo.
%% Now that corrections to old sartopo flybies have been implemented it may be better to average the two....
rx=2575+h2;
ry=2575+h2;
rz=2575+h2;

  

%%% Initial (t=0) pole right ascension, declination, rotation, and rotation rate
ra0=model.ra;   % units are deg
dec0=model.dec;  % units are deg
theta0=model.pm; % units are deg
w0=model.w; % units are deg/day





N=length(t2); % number of tiepoints



%% initialize resid structure
resid.wlon1=zeros(1,N);
resid.wlon2=zeros(1,N);
resid.lat1=zeros(1,N);
resid.lat2=zeros(1,N);
resid.dist=zeros(1,N);

% compute Euler angles, rotation matrices, and their derivatives and update matrix elements for each tiepoint
for i = 1:N
% Compute
% Euler angles rot mats and their derivatives for each observation
cent1=t1(i)/(365.25*100*24*3600);
cent2=t2(i)/(365.25*100*24*3600);

d1=t1(i)*second_to_day;
d2=t2(i)*second_to_day;
dd=1/24.0;
d1p=d1+dd; % used to compute instantaneous spin rate
d2p=d2+dd; % used to compute instantaneous spin rate



%% compute three eurler angles: alpha, beta, and theta for both obersvations of each tiepoint
alpha1=pi/2-dec0*dtr;
beta1=pi/2+ra0*dtr;
alpha2=pi/2-dec0*dtr;
beta2=pi/2+ra0*dtr;

%%% if there are long term pole trends in declination (ddec) and roight
%%% acsension (dra) include those in alphas and betas
if(isfield(model,'ddec'))

   alpha1=alpha1-model.ddec*cent1*dtr;
   beta1=beta1+model.dra*cent1*dtr;
   alpha2=alpha2-model.ddec*cent2*dtr;
   beta2=beta2+model.dra*cent2*dtr;
end  

%%% include constant spin rate term in thetas
theta1=theta0*dtr + w0*d1*dtr;
theta2=theta0*dtr + w0*d2*dtr;
theta1p=theta0*dtr + w0*d1p*dtr;
theta2p=theta0*dtr + w0*d2p*dtr;


%% incorporate periodic spin rate terms in thetas
for c=1:model.N_spin_period
c1w=cos(2*pi*d1/model.spin_period(c));
s1w=sin(2*pi*d1/model.spin_period(c));
c2w=cos(2*pi*d2/model.spin_period(c));
s2w=sin(2*pi*d2/model.spin_period(c));
c1wp=cos(2*pi*d1p/model.spin_period(c));
s1wp=sin(2*pi*d1p/model.spin_period(c));
c2wp=cos(2*pi*d2p/model.spin_period(c));
s2wp=sin(2*pi*d2p/model.spin_period(c));
theta1=theta1+dtr*model.spin_acos(c)*c1w+dtr*model.spin_asin(c)*s1w;
theta2=theta2+dtr*model.spin_acos(c)*c2w+dtr*model.spin_asin(c)*s2w;
theta1p=theta1p+dtr*model.spin_acos(c)*c1wp+dtr*model.spin_asin(c)*s1wp;
theta2p=theta2p+dtr*model.spin_acos(c)*c2wp+dtr*model.spin_asin(c)*s2wp;
end

%%% incorporate periodic pole declination and right ascension terms in
%%% alphas and betas
for c=1:model.N_pole_period
c1w=cos(2*pi*cent1/model.pole_period(c));
s1w=sin(2*pi*cent1/model.pole_period(c));
c2w=cos(2*pi*cent2/model.pole_period(c));
s2w=sin(2*pi*cent2/model.pole_period(c));
alpha1=alpha1-dtr*model.dec_acos(c)*c1w-dtr*model.dec_asin(c)*s1w;
alpha2=alpha2-dtr*model.dec_acos(c)*c2w-dtr*model.dec_asin(c)*s2w;
beta1=beta1+dtr*model.ra_acos(c)*c1w+dtr*model.ra_asin(c)*s1w;
beta2=beta2+dtr*model.ra_acos(c)*c2w+dtr*model.ra_asin(c)*s2w;
end



%rotate desired pole location to prime meridian
M3t1=[cos(beta1),sin(beta1),0;-sin(beta1),cos(beta1),0;0,0,1];
M3t2=[cos(beta2),sin(beta2),0;-sin(beta2),cos(beta2),0;0,0,1];

%rotate pole to z
M2t1=[1,0,0;0,cos(alpha1),sin(alpha1);0,-sin(alpha1),cos(alpha1)];
M2t2=[1,0,0;0,cos(alpha2),sin(alpha2);0,-sin(alpha2),cos(alpha2)];

%rotate prime meridian to its correct location
% this generates the inertial to titan body fixed rotation matrices
% Mt1 and Mt2 for times t1 and t2 respectively.
M1t1=[cos(theta1),sin(theta1),0;-sin(theta1),cos(theta1),0;0,0,1];
Mt1=M1t1*M2t1*M3t1;
M1t2=[cos(theta2),sin(theta2),0;-sin(theta2),cos(theta2),0;0,0,1];
Mt2=M1t2*M2t2*M3t2;




%% compute instanteous spin rates
w1=rtd*(theta1p-theta1)/dd;
w2=rtd*(theta2p-theta2)/dd;
     


%% compute inertial coordinates of tiepoints from raw measurement data and spin rate
%% spin rate only has a small second order effect
xi1(:,i)=get_j2000_pos(rx(i),ry(i),rz(i),Mt1,w1,dop1(i),range1(i),scpos1(:,i),scvel1(:,i),dc1(:,i),lambda1(i));
xi2(:,i)=get_j2000_pos(rx(i),ry(i),rz(i),Mt2,w2,dop2(i),range2(i),scpos2(:,i),scvel2(:,i),dc2(:,i),lambda2(i));

%% compute inertial coordinates if one assumes tiepoints are on the 2575-km radius reference sphere
xi1zeroh(:,i)=get_j2000_pos(2575,2575,2575,Mt1,w1,dop1(i),range1(i),scpos1(:,i),scvel1(:,i),dc1(:,i),lambda1(i));
xi2zeroh(:,i)=get_j2000_pos(2575,2575,2575,Mt2,w2,dop2(i),range2(i),scpos2(:,i),scvel2(:,i),dc2(:,i),lambda2(i));

%% compute body fixed coordinates using rotation matrices
xb1=M1t1*M2t1*M3t1*xi1(:,i);
xb2=M1t2*M2t2*M3t2*xi2(:,i);

%% compute body fixed coordinates if we assume tiepoints are on the reference sphere
xb1zeroh=M1t1*M2t1*M3t1*xi1zeroh(:,i);
xb2zeroh=M1t2*M2t2*M3t2*xi2zeroh(:,i);

%% compute residual position differences between observations of each tiepoint
dxb=xb1-xb2;
dxbzeroh=xb1zeroh-xb2zeroh;

%%% save intermediate values of first tiepoint for debugging
if(i==1)
 debug.dop1=dop1(1);
 debug.rx=rx(1);
 debug.ry=ry(1);
 debug.rz=rz(1);
 debug.range1=range1(1);
debug.scpos1=scpos1(:,1);
debug.scvel1=scvel1(:,1);
debug.dc1=dc1(:,1);
 debug.lambda1=lambda1(1);


 debug.dop2=dop2(1);
 debug.range2=range2(1);
debug.scpos2=scpos2(:,1);
debug.scvel2=scvel2(:,1);
debug.dc2=dc2(:,1);
 debug.lambda2=lambda2(1);

debug.t1=t1(1);
debug.t2=t2(1);
  debug.alpha1=alpha1;
debug.beta1=beta1;
debug.theta1=theta1;
  debug.alpha2=alpha2;
debug.beta2=beta2;
debug.theta2=theta2;

debug.Mt1=Mt1;
debug.Mt2=Mt2;
debug.xi1=xi1(:,1);
debug.xi2=xi2(:,1);
debug.xb1=xb1;
debug.xb2=xb2;
debug.w1=w1;
debug.w2=w2;
end

%% compute residual distances and body fixed longitudes and latitudes
%% for SARTopo height and reference sphere (height=0) cases
%% for both observations
resid.dist(i)=sqrt(sum(dxb.*dxb));
distzeroh(i)=sqrt(sum(dxbzeroh.*dxbzeroh));
[resid.wlon1(i),resid.lat1(i)]=pos_to_lonlat(xb1);
[resid.wlon2(i),resid.lat2(i)]=pos_to_lonlat(xb2);

[lon1zeroh(i),lat1zeroh(i)]=pos_to_lonlat(xb1zeroh);
[lon2zeroh(i),lat2zeroh(i)]=pos_to_lonlat(xb2zeroh);
end

resid.rms_dist=rms(resid.dist);

%%% Put together output tiepts structures
tiepts.name=names(gg);
tiepts.height=h2;

tiepts.time_1=t1;
tiepts.doppler_1=dop1;
tiepts.range_1=range1;
tiepts.scpos_inert_1=scpos1';
tiepts.scvel_inert_1=scvel1';
tiepts.lookvec_inert_1=(2575.0)*dc1'-scpos1'; % dc1 is a unit vector
tiepts.wavelength_1=lambda1;
tiepts.pos_inert_1=xi1';
tiepts.lat_1=resid.lat1;
tiepts.wlon_1=resid.wlon1;
tiepts.pos_inert_zeroh_1=xi1zeroh';
tiepts.lat_zeroh_1=lat1zeroh;
tiepts.wlon_zeroh_1=lon1zeroh;


tiepts.time_2=t2;
tiepts.doppler_2=dop2;
tiepts.range_2=range2;
tiepts.scpos_inert_2=scpos2';
tiepts.scvel_inert_2=scvel2';
tiepts.lookvec_inert_2=(2575.0)*dc2'-scpos2';
tiepts.wavelength_2=lambda2;
tiepts.pos_inert_2=xi2';
tiepts.lat_2=resid.lat2;
tiepts.wlon_2=resid.wlon2;
tiepts.pos_inert_zeroh_2=xi2zeroh';
tiepts.lat_zeroh_2=lat2zeroh;
tiepts.wlon_zeroh_2=lon2zeroh;
tiepts.residual_dist=resid.dist;
tiepts.residual_dist_zeroh=distzeroh;

