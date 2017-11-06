function out = MTLB_name(name,num)
out = [name,regexprep(num2str(num),'\s+','_')];

end