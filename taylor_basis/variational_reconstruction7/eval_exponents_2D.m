function [exp,val,val_compare] = eval_exponents_2D(degree,pt)
N = nchoosek(2+degree,degree);
exp = zeros(2,N);
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

for j = 0:degree
    for i = 0:degree-j
        cnt = cnt + 1;
        exp(1,cnt) = i;
        exp(2,cnt) = j;
        for ii = 1:i
            val(cnt) = val(cnt)*pt(1);
        end
        for jj = 1:j
            val(cnt) = val(cnt)*pt(2);
        end
    end
end

val_compare = pt(1).^exp(1,:) .* pt(2).^exp(2,:);

end