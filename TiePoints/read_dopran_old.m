function [dop1,range1,scpos1,scvel1,dc1,t1,lambda1,dop2,range2,scpos2,scvel2,dc2,t2,lambda2,h2,names,sab1,sab2]=read_dopran (filename);
%%% routine to read values from a raw tiepoints file

fid=fopen(filename,'r');
i=1;
while ~feof(fid) 
  str=fscanf(fid,'%s',[1,1]);
  names{i}=str;
  if(~feof(fid))
  range1(:,i)=fscanf(fid,'%g',[1,1]);
  dop1(:,i)=fscanf(fid,'%g',[1,1]);
  scpos1(:,i)=fscanf(fid,'%g',[3,1]);
  scvel1(:,i)=fscanf(fid,'%g',[3,1]);
  dc1(:,i)=fscanf(fid,'%g',[3,1]);
  t1(:,i)=fscanf(fid,'%g',[1,1]);  
  lambda1(:,i)=fscanf(fid,'%g',[1,1]);  
  range2(:,i)=fscanf(fid,'%g',[1,1]);
  dop2(:,i)=fscanf(fid,'%g',[1,1]);
  scpos2(:,i)=fscanf(fid,'%g',[3,1]);
  scvel2(:,i)=fscanf(fid,'%g',[3,1]);
  dc2(:,i)=fscanf(fid,'%g',[3,1]);
  t2(:,i)=fscanf(fid,'%g',[1,1]);  
  lambda2(:,i)=fscanf(fid,'%g',[1,1]);  
  h2(:,i)=fscanf(fid,'%g',[1,1]);  
  %sab1(:,i)=fscanf(fid,'%d',[1,1]);  
  %sab2(:,i)=fscanf(fid,'%d',[1,1]);  
  end
  i=i+1;
end
fclose(fid);
