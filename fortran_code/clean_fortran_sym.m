function str = clean_fortran_sym(expr)
prec_str = '_dp';
max_chars = 75;
constant_dict = {'3.141592653589793D0', 'pi', ...
                 '0.0D0', 'zero' , ...
                 '1.0D0', 'one'  , ...
                 '2.0D0', 'two'  , ...
                 '3.0D0', 'three', ...
                 '4.0D0', 'four' , ...
                 '5.0D0', 'five' , ...
                 '6.0D0', 'six'  , ...
                 '7.0D0', 'seven', ...
                 '8.0D0', 'eight', ...
                 '9.0D0', 'nine' , ...
                 '1.0D+1', 'ten' };
constant_dict = reshape(constant_dict,2,[]);
str = fortran(expr);

str = strip(str,'left');
str = strrep(str,newline,'');
str = regexprep(str, '\s*&', '');
str = replace(str,constant_dict(1,:),constant_dict(2,:));
% str = regexprep(str, '(?<=(\d+\.\d*|\.\d+|\d+))D(?=[+-]?\d+)', 'e');
str = regexprep(str, '(\d+\.\d*|\.\d+|\d+)D([+-]?\d+)', ['$1e$2',prec_str]);
str = regexprep(str, '(\d+)_(\d+)', '($1,$2)');
str = regexprep(str, '([a-zA-Z])(\d+)', '$1($2)');
% str = regexprep(str, '(?<!\*)\*(?!\*)', ' $& ');
str = regexprep(str, '[+=]', ' $& ');
str = regexprep(str, '(?<!\*)\*(?!\*)|(?<!-)(?<![eE])-(?!-)|(?<!-)-(?!-)(?!\d)', ' $& ');

wrappedCell = textwrap({str}, max_chars);

str = sprintf('%s &\n', wrappedCell{:});

str = regexprep(str, '&\s*$', '');
str = strip(str, 'right');
end