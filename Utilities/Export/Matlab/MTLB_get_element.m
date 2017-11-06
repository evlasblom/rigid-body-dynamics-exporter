function out = MTLB_get_element(name,num,in)
% out = {name,numbers,caller,expression}
out.name = name;
out.num = in.num;
out.caller = MTLB_name(name,[num{1},0,num{2}]);

% select elements (num should be a 1x2 cell)
zero = false;
str = [in.caller,'('];
for ii = 1:2
    str = [str,'[',num2str(num{ii}),'],'];
end
str = [str(1:end-1),')'];
out.exp = str;

% zero
out.zero = zero;

end