% Converts symbolic matrix into an anonymous function.
% Used in Animation3D.m
function f = HomogeneousTransformFunction(H)

if isa(H,'sym') % is symbolic? get inputs, use char
    qlocal = symvar(H);
    fnc_input = '';
    if ~isempty(qlocal)
        for ii = 1:length(qlocal)
            fnc_input = [fnc_input,char(qlocal(ii)),','];
        end
    end
    fnc_convert = @char;
else % not symbolic? no inputs
    fnc_input = '';
    fnc_convert = @num2str;
end

string = ['@(',fnc_input(1:end-1),') '];
for ii = 1:size(H,1)
    if ii == 1; string=[string,'[']; end
    for jj = 1:size(H,2)
        string = [string,fnc_convert(H(ii,jj)),','];
        if jj == size(H,2); string = [string(1:end-1),';']; end
    end
    if ii == size(H,1); string = [string(1:end-1),'];']; end
end

f = str2func(string);



