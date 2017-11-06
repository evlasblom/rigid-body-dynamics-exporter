function out = MTLB_addition(name,num,varargin)
% out = {name,numbers,caller,expression}
out.name = name;
out.num = num;
out.caller = MTLB_name(name,num);

zero = true;
str = '(';
for ii = 1:nargin-2
    if varargin{ii}.zero
        continue
    else
        zero = false;
        if isempty(varargin{ii}.caller); term = varargin{ii}.exp;
        else term = varargin{ii}.caller;
        end
        str = [str,term,'+'];
    end
end
str(end) = ')';
out.exp = str;

% zero
out.zero = zero;

end