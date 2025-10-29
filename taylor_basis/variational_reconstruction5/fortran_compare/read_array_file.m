function A = read_array_file( file_name )

fid = fopen(file_name,'r');


line = fgetl(fid);
dim = regexp(line,'\d*','match');
dim = cellfun(@str2num,dim);
A = fscanf(fid,'%f');

A = reshape(A,[dim,1]);
fclose(fid);

end