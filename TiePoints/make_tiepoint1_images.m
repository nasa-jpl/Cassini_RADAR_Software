% produce images of the first tiepoint in the good files

addpath ~/cassini/matlab

id='tas01t23s01_jpl2_p1';

lat1=51.1331;
wlon1=80.5191;
lat2=51.0332;
wlon2=80.4920;
fs0_1='BIBQI49N071_D035_T00AS01_V03.IMG';
fs0_2='BIBQI26N009_D111_T023S01_V03.IMG';
flat_1='BITQI49N071_D035_T00AS01_V03.IMG';
flat_2='BITQI26N009_D111_T023S01_V03.IMG';
fwlon_1='BINQI49N071_D035_T00AS01_V03.IMG';
fwlon_2='BINQI26N009_D111_T023S01_V03.IMG';

reload=0;
if(reload)
lines=26368;
samples=4096;
pixtype='uchar';
n_headrec=2;
s0arr1=read_bidr(fs0_1,lines,samples,pixtype,n_headrec);

pixtype='float';
n_headrec=1;
latarr1=read_bidr(flat_1,lines,samples,pixtype,n_headrec);
wlonarr1=read_bidr(fwlon_1,lines,samples,pixtype,n_headrec);


lines=38400;
samples=3840;
pixtype='uchar';
n_headrec=2;
s0arr2=read_bidr(fs0_2,lines,samples,pixtype,n_headrec);

pixtype='float';
n_headrec=1;
latarr2=read_bidr(flat_2,lines,samples,pixtype,n_headrec);
wlonarr2=read_bidr(fwlon_2,lines,samples,pixtype,n_headrec);

dlat1=latarr1-lat1;
dlon1=(wlonarr1-wlon1)*cos(lat1*pi/180);
dlat2=latarr2-lat2;
dlon2=(wlonarr2-wlon2)*cos(lat2*pi/180);
dist1=sqrt(dlat1.*dlat1+dlon1.*dlon1);
dist2=sqrt(dlat2.*dlat2+dlon2.*dlon2);
[mdist1,idx1]=min(dist1(:));
[mdist2,idx2]=min(dist2(:));
sz1=size(s0arr1);
sz2=size(s0arr2);
[i1,j1]=ind2sub(sz1,idx1);
[i2,j2]=ind2sub(sz2,idx2);

end

rad=250;
im1=s0arr1((i1-rad):(i1+rad),(j1-rad):(j1+rad));
im2=s0arr2((i2-rad):(i2+rad),(j2-rad):(j2+rad));

figure(1)
colormap('gray');
y=(i1-rad):(i1+rad);
x=(j1-rad):(j1+rad);
imagesc(x,y,im1);colorbar;
hold on
h=plot(j1,i1,'rx');
set(h,'MarkerSize',50);
set(h,'LineWidth',2);
hold off
xlabel(sprintf('Sample number (Nearest=%d)',j1));
ylabel(sprintf('Line number (Nearest=%d)',i1));
title(sprintf('Tiepoint %s as observed in TA flyby',id),'Interpreter','none');

print('-dpng',sprintf('TiepointImage_%s_TA.png',id));


figure(2)
colormap('gray');
y=(i2-rad):(i2+rad);
x=(j2-rad):(j2+rad);
imagesc(x,y,im2);colorbar;
hold on
h=plot(j2,i2,'rx');
set(h,'MarkerSize',50);
set(h,'LineWidth',2);
hold off
xlabel(sprintf('Sample number (Nearest=%d)',j2));
ylabel(sprintf('Line number (Nearest=%d)',i2));
title(sprintf('Tiepoint %s as observed in T23 flyby',id),'Interpreter','none');

print('-dpng',sprintf('TiepointImage_%s_T23.png',id));
