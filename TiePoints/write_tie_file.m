function write_tie_file(fname,tiepts,outfields,fmt)
%%function write_tie_file(fname,tiepts,outfields,fmt)
%% function to write tiepoint file given
%% fname = filename of output file
%% tiepts = tiepoint information structure
%% outfields = list of fields from structure to write
%% fmt = format to use in fprintf call

%% open output file
  ofp=fopen(fname,'w');

%% get lengths of filename and numebr of fields to output  
  N=length(tiepts.name);
  M=length(outfields);

%% looping over list of fields, assign values to output
%% and convert field names to strings for column headings
  vals=[];
  off=1;
  for j=1:M
     v=getfield(tiepts,outfields{j});
     sz=prod(size(v));
     if(sz==N)
       v=reshape(v,[N,1]);
       outlabels{off}=outfields{j};
       off=off+1;
     else
       v=reshape(v,[N,3]);
       outlabels{off}=sprintf('%s_x',outfields{j});
       off=off+1;
       outlabels{off}=sprintf('%s_y',outfields{j});
       off=off+1;
       outlabels{off}=sprintf('%s_z',outfields{j});
       off=off+1;
     end
     vals=cat(2,vals,v);
  end

  %%% Write column headers to file
  fprintf(ofp,'Tiepoint_ID ');
  fprintf(ofp,',%s',outlabels{:});
  fprintf(ofp,'\n');

  %%% loop over tiepoints and write each row to output file
  for i=1:N
    fprintf(ofp,fmt,tiepts.name{i},vals(i,:));
  end

  %%% close file
  fclose(ofp);

