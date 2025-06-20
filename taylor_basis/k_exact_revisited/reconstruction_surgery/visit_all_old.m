function [val,list] = visit_all_old(list,order)
N = numel(list);
val = list(order(N));
list(order(N)) = 0;
if val == 0
    for i = N-1:1
        if list(order(i)) ~= 0
            val = list(order(i));
            list(order(i)) = 0;
            return
        end
    end
    error('value not in list');
end

end