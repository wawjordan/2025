function IDX = generate_index_array(idx_start,idx_end)
idx_hi = idx_end - idx_start + 1;
tmp = num2cell(1:prod(idx_hi));
lin_idx = reshape(sub2ind(idx_hi,[tmp{:}]),idx_hi);
IDX = arrayfun(@(n)my_ind2sub(idx_hi,n) + idx_start - 1,lin_idx,UniformOutput=false);
end

function idx_out = my_ind2sub(sz,ndx)
% [v1,v2,varargout] = ind2sub(siz,ndx)
[idx_tmp{1:numel(sz)}] = ind2sub(sz,ndx);
idx_out = cell2mat(idx_tmp);
end