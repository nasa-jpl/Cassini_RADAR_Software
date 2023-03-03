function im = read_bidr(file,lines,samples, pixtype, num_header_recs)

% for now takes lines and samples , eventually should reader header
fid=fopen(file,'r','l');

% read past header commented out for now
dummy=fread(fid,[samples,num_header_recs],pixtype);

im=fread(fid,[samples,lines],pixtype);
im=im';

