% m = 2; cnt = 0; for p = 0:m, fprintf('total order=%d\n',p); for k=0:p, for l=0:k, cnt = cnt+1; fprintf('%-.2d: %d + %d + %d = %d\n',cnt,p-k,k-l,l,p); end, end, end

function exp = my_grevlex(m,dim)
dim = 3;
exp = zeros(nchoosek(m+dim,dim),dim);
cnt = 0;
% for p = 0:m
%     fprintf('total order=%d\n',p);
%     for k=0:p
%         for l=0:k
%             cnt = cnt+1;
%             exp(cnt,:) = [p-k,k-l,l];
%             exp_list = num2cell(exp(cnt,:));
%             fprintf(['%-.2d: ',repmat('%d + ',[1,dim-1]),'%d = %d\n'],cnt,exp_list{:},p);
%         end
%     end
% end

exp = exp_loop(m,dim);

end

function exp = exp_loop(m,dim)
n_terms = nchoosek(m+dim,dim);
exp = zeros(n_terms,dim);
cnt = 0;
% base level
for nd = 2:dim
    for j = 0:m
        idx = zeros(1,dim);
        idx(nd) = j;
        for d = 2:dim
            for i = 0:j-(d-1)
                idx(d) = i;
                cnt = cnt + 1;
                fprintf('%d:%d,%d\n',cnt,d,j);
                % fprintf('%d:%d,%d\n',cnt,d,j);
            end
        end
    end
end
end

% function exp = exp_loop(m,dim)
% n_terms = nchoosek(m+dim,dim);
% exp = zeros(n_terms,dim);
% cnt = 0;
% % base level
% for j = 0:m
%     idx = zeros(1,dim);
%     for d = 1:dim
%         for i = 0:j-(d-1)
%             cnt = cnt + 1;
%             idx(d) = j+i;
%             fprintf('%d:%d,%d,%d\n',cnt,idx);
%             % fprintf('%d:%d,%d\n',cnt,d,i);
%         end
%     end
% end
% end

function iG = local_to_global(iSub,nSub)
iSub = iSub(:);
nSub = nSub(:);
nDims = numel(iSub);
p = 1;
iG = 1;
for i = 1:nDims
    iG = iG + ( iSub(i) - 1 )*p;
    p = p*nSub(i);
end
end

function iSub = global_to_local(iG,nSub)
nDims = numel(nSub);
iSub = zeros(1,nDims);
if (nDims==1)
    iSub(1) = iG;
    return
end
p = prod(nSub);
iGtmp = iG;
for i = nDims:-1:1
    p = fix( p/nSub(i) );
    iTmp = mod(iGtmp-1,p)+1;
    iSub(i) = fix( (iGtmp-iTmp)/p ) + 1;
    iGtmp = iTmp;
end
end