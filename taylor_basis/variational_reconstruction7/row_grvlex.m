function [e_out,i] = row_grvlex(e_in)
[~,i] = sortrows([sum(e_in,2),fliplr(e_in)]);
e_out = e_in(i,:);
end