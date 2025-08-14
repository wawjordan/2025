function [exp,val,val_compare] = eval_exponents_3D(degree,pt)
N = nchoosek(3+degree,degree);
exp = zeros(3,N);
val = ones(1,N);
cnt = 0;
% for j = 0:degree
%     for i = j:degree
%         cnt = cnt + 1;
%         exp(1,cnt) = j;
%         exp(2,cnt) = i-j;
%         for jj = 1:j
%             val(cnt) = val(cnt)*pt(1);
%         end
%         for ii = 1:i-j
%             val(cnt) = val(cnt)*pt(2);
%         end
%     end
% end

for k = 0:degree
    for j = 0:degree-k
        for i = 0:degree-j-k
            cnt = cnt + 1;
            fprintf('%d %d %d %d\n',cnt,i,j,k)
            exp(1,cnt) = i;
            exp(2,cnt) = j;
            exp(3,cnt) = k;
            for ii = 1:i
                val(cnt) = val(cnt)*pt(1);
            end
            for jj = 1:j
                val(cnt) = val(cnt)*pt(2);
            end
            for kk = 1:k
                val(cnt) = val(cnt)*pt(3);
            end
        end
    end
end
val_compare = ( pt(1).^exp(1,:) ) .* ( pt(2).^exp(2,:) ) .* ( pt(3).^exp(3,:) );

end