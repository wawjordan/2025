function [n_faces,face_indices] = boundary_stencil_intersect_old(dim,stencil,grid)
stencil_size = size(stencil,2);
% bnd_min = zeros(4,grid.nblocks);
% bnd_max = zeros(4,grid.nblocks);
% in_block = false(grid.nblocks,1);

max_bcs = 2*dim;

% for n = 1:grid.nblocks
%     max_bcs = max(max_bcs,grid.gblock(n).nbounds);
%     [bnd_min(:,n),bnd_max(:,n),in_block(n)] = get_bounding_box(stencil,n);
% end

% face_mask = false(6,stencil_size);
face_mask = false(max_bcs,stencil_size);
for s = 1:stencil_size
    block = stencil(1,s);
    lo_bnd_g = [1;1;1];
    hi_bnd_g = grid.gblock(block).Ncells;
    idx = stencil(2:4,s);
    % for n = 1:grid.gblock(block).nbounds
    % for n = 1:6
    %     face_mask(n,s) = in_bound(idx,lo_bnd_g,hi_bnd_g);
    % end
    face_mask(1:2*dim,s) = on_nd_boundary(dim,idx,lo_bnd_g,hi_bnd_g);
end
n_faces = sum(face_mask,"all");

face_indices = zeros(4,n_faces);
cnt = 0;
for s = 1:stencil_size
    % for n = 1:grid.gblock(block).nbounds
    for n = 1:2*dim
        if face_mask(n,s)
            cnt = cnt + 1;
            face_indices(1,cnt)   = stencil(1,s);
            face_indices(2:4,cnt) = stencil(2:4,s);
            tmp_idx = fix((n+1)/2)+1;
            face_indices(tmp_idx,cnt) = face_indices(tmp_idx,cnt) + mod(n+1,2);
        end
    end
end

end

% function val = cell_offsets()
% val = zeros([3,6],"int32");
% direction = [-1,1];
% for d = 1:3
%     for s = 1:2
%         val(d,2*(d-1)+s) = direction(s);
%     end
% end
% end
% for i = 1:6, d = [0,0,0]; d(fix((i+1)/2)) = mod(i+1,2); fprintf('%d:[%d,%d]->[%d,%d,%d]\n',i,fix((i+1)/2),mod(i+1,2),d); end

% function [face_indices,n_faces] = boundary_stencil_intersect_old(stencil,mask,grid)
% stencil_size = size(stencil,2);
% bnd_min = zeros(4,grid.nblocks);
% bnd_max = zeros(4,grid.nblocks);
% in_block = false(grid.nblocks,1);
% 
% max_bcs = 0;
% 
% for n = 1:grid.nblocks
%     max_bcs = max(max_bcs,grid.gblock(n).nbounds);
%     [bnd_min(:,n),bnd_max(:,n),in_block(n)] = get_bounding_box(stencil,n);
% end
% 
% % face_mask = false(6,stencil_size);
% face_mask = false(max_bcs,stencil_size);
% for s = 1:stencil_size
%     if mask(s)
%         continue
%     end
%     block = stencil(1,s);
%     if ~in_block(block)
%         continue
%     end
%     lo_bnd_s = bnd_min(2:4,block);
%     hi_bnd_s = bnd_max(2:4,block);
%     lo_bnd_g = [1;1;1];
%     hi_bnd_g = grid.gblock(block).Ncells;
%     if ~bound_intersect( lo_bnd_g, hi_bnd_g, lo_bnd_s, hi_bnd_s )
%         continue
%     end
%     idx = stencil(2:4,s);
%     % face_mask(:,s) = on_3d_boundary(idx,lo_bnd_g,hi_bnd_g);
%     for n = 1:grid.gblock(block).nbounds
%         face_mask(n,s) = in_bound(idx,lo_bnd_g,hi_bnd_g);
%     end
% end
% n_faces = sum(face_mask,[],"all");
% 
% face_indices = zeros(4,n_faces);
% cnt = 0;
% for s = 1:stencil_size
%     for n = 1:grid.gblock(block).nbounds
%         if face_mask(n,s)
%             cnt = cnt + 1;
%             face_indices(1,cnt)   = stencil(1,s);
%             face_indices(2:4,cnt) = stencil(2:4,s) + grid.gblock(block).bound.out_dir(:);
%         end
%     end
% end
% 
% end

% function [bnd_min,bnd_max,in_block] = get_bounding_box(stencil,block)
% bnd_min = zeros(4,1); bnd_min(1) = block;
% bnd_max = zeros(4,1); bnd_max(1) = block;
% in_block = false;
% mask = stencil(1,:)==block;
% n_mask = sum(mask);
% if n_mask == 0
%     return
% end
% 
% bnd_min = min(stencil(:,mask),[],2);
% bnd_max = max(stencil(:,mask),[],2);
% end
% 
% function val = bound_intersect(bnd1_min, bnd1_max, bnd2_min, bnd2_max)
% val = all( range_intersect( bnd1_min, bnd1_max, bnd2_min, bnd2_max ) );
% end
% 
% function val = range_intersect( startA, endA, startB, endB )
%   val = ( startA <= endB )&&( endA >= startB );
% end
% 
% function val = on_3d_boundary(idx,bnd_min,bnd_max)
% val = false(6,1);
% bnds = [bnd_min(:),bnd_max(:)];
% for d = 1:3
%     for s = 1:2
%         bnds_tmp = bnds;
%         bnds_tmp(d,:) = bnds(d,s);
%         val(2*(d-1)+s) = in_bound(idx,bnds_tmp(:,1),bnds_tmp(:,2));
%     end
% end
% end

function val = on_nd_boundary(n_dim,idx,bnd_min,bnd_max)
val = false(n_dim*2,1);
bnds = [bnd_min(:),bnd_max(:)];
for d = 1:n_dim
    for s = 1:2
        bnds_tmp = bnds;
        bnds_tmp(d,:) = bnds(d,s);
        val(2*(d-1)+s) = in_bound(idx,bnds_tmp(:,1),bnds_tmp(:,2));
    end
end
end

function val = in_bound(idx,bnd_min,bnd_max)
    val1 = all(idx(:)>=bnd_min(:))&&all(idx(:)<=bnd_max(:));
    val2 = all(idx(:)<=bnd_min(:))&&all(idx(:)>=bnd_max(:));
    val = val1||val2;
end
