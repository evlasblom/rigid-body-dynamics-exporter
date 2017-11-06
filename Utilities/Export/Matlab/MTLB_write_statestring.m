function out = MTLB_write_statestring(allstates)

str = '';
for i = 1:length(allstates)
    str = [str,char(allstates(i)),' = q(',num2str(i),'); '];
end
str = str(1:end-2);
str = [str,';'];

out = sprintf('%s\n\n',str);

end