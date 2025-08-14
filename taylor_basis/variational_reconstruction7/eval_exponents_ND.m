function [N] = eval_exponents_ND(dim,degree)
% function [exp,val,val_compare] = eval_exponents_ND(dim,degree,pt)
N = nchoosek(dim+degree,degree);
exp = zeros(dim,N);
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

num_fmt = '%5d';
fmt_fun = @(n,fmt) [repmat([fmt,' '],1,n-1), fmt, '\n'];

fmt = fmt_fun(dim+1,num_fmt);
% fmt = fmt_fun(1,num_fmt);

sum_list = zeros(dim,degree);
for d = 1:dim
    for n = 0:degree
    sum_list()

for d = 1:dim
    exp0(d) = exp0(d) + 1;
    exp_sum = sum(exp0);
    for i = 0:degree-exp_sum

        for in = 0:
            cnt = cnt + 1;
            fprintf(fmt,[dims;cnt]);
        end
    end
end

% fprintf(fmt,(,cnt,i,j,k)

end