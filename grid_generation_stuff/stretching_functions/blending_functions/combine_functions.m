function f = combine_functions(f,w)
f = @(t) arrayfun(@(t)combined_function_value(t,f,w),t);
end
function val = combined_function_value(t,f,w)
    n = numel(f);
    den = w{1}(t);
    for i = 2:n
        den = den + w{i}(t);
    end
    val = f{1}(t).*w{1}(t)./den;
    for i = 2:n
        val = val + f{i}(t).*w{i}(t)./den;
    end
end