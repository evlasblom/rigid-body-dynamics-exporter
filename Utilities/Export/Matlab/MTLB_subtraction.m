function out = MTLB_subtraction(name,num,varargin)
% out = {name,numbers,caller,expression}
out.name = name;
out.num = num;
out.caller = MTLB_name(name,num);

addminus = false;
zero = true;
str = '(';
for ii = 1:nargin-2
    if varargin{ii}.zero
        if ii == 1
            addminus = true;
        end
        continue
    else
        zero = false;
        if isempty(varargin{ii}.caller); term = varargin{ii}.exp;
        else term = varargin{ii}.caller;
        end
        str = [str,term,'-'];
    end
end
if addminus
    str = ['(','-',str(2:end)];
end
str(end) = ')';
out.exp = str;

% zero
out.zero = zero;

end