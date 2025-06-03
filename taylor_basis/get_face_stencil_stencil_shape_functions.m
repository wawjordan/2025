function [TS,QS] = get_face_stencil_stencil_shape_functions(GRID,blk,idx,degree)

TS = struct();
QS = struct();
[TS(1).T,QS(1).Q,idxs] = get_face_stencil_shape_functions(GRID,blk,idx,degree);
n_nbors = size(idxs,2);
for i = 2:n_nbors
    [TS(i).T,QS(i).Q,~] = get_face_stencil_shape_functions(GRID,blk,idxs(:,i),degree);
end

end