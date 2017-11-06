function out = MTLB_multiplication(name,num,varargin)
% out = {name,numbers,caller,expression}
out.name = name;
out.num = num;
out.caller = MTLB_name(name,num);

zero = false;
str = '';
for ii = 1:nargin-2
    if varargin{ii}.zero
        zero = true;
        break
    else
        if isempty(varargin{ii}.caller); term = varargin{ii}.exp;
        else term = varargin{ii}.caller;
        end
        str = [str,term,'*'];
    end
end
out.exp = str(1:end-1);

% zero
out.zero = zero;

end