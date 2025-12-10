function SENSEI_BC_write(N,idx_list,bc_id_list,additional_inputs,filename)

N_bcs = length(bc_id_list);

fid = fopen(filename,'w');
intfmt = ' %-5d';
% fltfmt = ' %-# 23.16E';
fltfmt = ' %-.16g';

fprintf(fid,'%d\n',1);
fprintf(fid,'%d\n',1);
fprintf(fid,[repmat(intfmt,[1,length(N)]),'\n'],N(:));
fprintf(fid,'A\n');
fprintf(fid,'%d\n',N_bcs);

for i = 1:N_bcs
    bc_name = sprintf('bc%d  ',i);
    fmt1 = repmat(intfmt,[1,length(idx_list{i})+1]);
    fmt2 = repmat(fltfmt,[1,length(additional_inputs{i})]);
    tmp = [idx_list{i}(:); bc_id_list(i); additional_inputs{i}(:)];
    fprintf( fid, [bc_name,fmt1,fmt2,'\n'], tmp );
end

fclose(fid);
end