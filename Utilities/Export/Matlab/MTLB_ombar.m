% operator
function out = MTLB_ombar(in)
% out = {name,numbers,caller,expression}
out.name = [];
out.num = [];
out.caller = MTLB_name(out.name,out.num);

% write expression
if isempty(in.caller); term = in.exp;
else term = in.caller;
end
exp = sprintf('%s',['f_ombar(',term,')']);
out.exp = exp;

% zero
out.zero = in.zero;

end