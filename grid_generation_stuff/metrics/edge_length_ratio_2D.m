function L = edge_length_ratio_2D(x,y,dir)
% ratio of distance between consecutive points on a coordinate line
%  j-1     j     j+1
%   +------+------+
%   |<-L1->|<-L2->|
% L = min(L1,L2)/max(L1,L2)
% output is same size as input coordinates
% endpoint ratios are set to 1
sz = size(x);
imax = sz(1);
jmax = sz(2);

L = ones(sz);

switch dir
    case(1)
        L_tmp = sqrt( diff(x,[],1).^2 + diff(y,[],1).^2 );
        L(2:imax-1,:) = min(L_tmp(1:imax-2,:),L_tmp(2:imax-1,:)) ...
                     ./ max(L_tmp(1:imax-2,:),L_tmp(2:imax-1,:));
    case(2)
        L_tmp = sqrt( diff(x,[],2).^2 + diff(y,[],2).^2 );
        L(:,2:jmax-1) = min(L_tmp(:,1:jmax-2),L_tmp(:,2:jmax-1)) ...
                     ./ max(L_tmp(:,1:jmax-2),L_tmp(:,2:jmax-1));
    otherwise
        error('dir can only be 1 (i-coordinates), or 2 (j-coordinates)')
end

end