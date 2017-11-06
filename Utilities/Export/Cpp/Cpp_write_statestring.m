function Cpp_write_statestring(fid,allstates)

str = 'double ';
for i = 1:length(allstates)
    str = [str,char(allstates(i)),' = q[',num2str(i),']; '];
end
str = str(1:end-2);
str = [str,';'];

fprintf(fid,'%s\n\n',str);

end