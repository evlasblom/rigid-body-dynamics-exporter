function out = MTLB_write_parameterstring(parameters)

pnames = fieldnames(parameters);

str = '';
for i = 1:length(pnames)
    str = [str,pnames{i},' = p.',pnames{i},'; '];
end
str = str(1:end-2);
str = [str,';'];

out = sprintf('%s\n\n',str);

end