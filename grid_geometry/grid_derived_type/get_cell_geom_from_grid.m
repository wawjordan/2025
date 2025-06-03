function [center,quad,nodes,face_quads,normals,nbor_center] = get_cell_geom_from_grid(GRID,blk,idx,dim)

bnd_min = [1;1;1];
bnd_max = [ GRID.gblock(blk).imax - 1;...
            GRID.gblock(blk).jmax - 1;...
            GRID.gblock(blk).kmax - 1 ];

center    = GRID.gblock(blk).grid_vars.cell_c(:,idx(1),idx(2),idx(3));
quad      = GRID.gblock(blk).grid_vars.quad(idx(1),idx(2),idx(3));
node_var={'x','y','z'};
quad_var={'xi_face_quad','eta_face_quad','zeta_face_quad'};
norm_var={'xi_n','eta_n','zeta_n'};
if dim == 2
    nodes      = zeros(3,4);
    face_quads = cell(1,4);
    normals    = cell(1,4);
    nbor_center = cell(1,4);
    off   = [0,1,1,0; ...
             0,0,1,1; ...
             0,0,0,0];
    for i = 1:3
        for j = 1:4
            nodes(i,j) = GRID.gblock(blk).(node_var{i})( idx(1) + off(1,j),...
                                                         idx(2) + off(2,j),...
                                                         idx(3) + off(3,j) );
        end
    end
    for i = 1:4
        face_idx =  fix((i-1)/2)+1;      % [1,1,2,2]
        offset   =  mod(i+1,2);          % [0,1]
        factor   =  2*offset-1;          % [-1, 1]
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+offset;
        face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}.v = normals{i}.v * factor;
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+factor;
        if in_bounds(tmp_idx,bnd_min,bnd_max)
            nbor_center{i}  = GRID.gblock(blk).grid_vars.cell_c(:,tmp_idx(1),tmp_idx(2),tmp_idx(3));
        end
    end
else
    nodes      = zeros(3,8);
    face_quads = cell(1,6);
    normals    = cell(1,6);
    nbor_center = cell(1,6);
    off   = [0,1,1,0,0,1,1,0; ...
             0,0,1,1,0,0,1,1; ...
             0,0,0,0,1,1,1,1];
    for i = 1:3
        for j = 1:8
            nodes(i,j) = GRID.gblock(blk).(node_var{i})( idx(1) + off(1,j),...
                                                         idx(2) + off(2,j),...
                                                         idx(3) + off(3,j) );
        end
    end
    for i = 1:6
        face_idx =  fix((i-1)/2)+1; % [1,1,2,2,3,3]
        offset   =  mod(i+1,2);     % [0,1]
        factor   =  2*offset-1;     % [-1, 1]
        tmp_idx  =  idx;
        tmp_idx(face_idx) = tmp_idx(face_idx)+offset;
        face_quads{i} = GRID.gblock(blk).grid_vars.(quad_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}    = GRID.gblock(blk).grid_vars.(norm_var{face_idx}) ...
                                                   ( tmp_idx(1),...
                                                     tmp_idx(2),...
                                                     tmp_idx(3) );
        normals{i}.v = normals{i}.v * factor;
        if in_bounds(tmp_idx,bnd_min,bnd_max)
            nbor_center{i}  = GRID.gblock(blk).grid_vars.cell_c(:,tmp_idx(1),tmp_idx(2),tmp_idx(3));
        end
    end
end
end

function val = in_bounds(idx,bnd_min,bnd_max)
val = ( all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:)) || ...
        all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:)) );
end