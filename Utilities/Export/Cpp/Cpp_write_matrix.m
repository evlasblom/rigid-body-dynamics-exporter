function Cpp_write_matrix(fid,Mname,M)

s1 = size(M,1);
s2 = size(M,2);

% gsl_matrix M
fprintf(fid,'%s',...
    ['gsl_matrix* ',Mname,' = gsl_matrix_alloc(',num2str(s1),num2str(s2),');']);


% how to do this?

% [a,b,c; d,e,f; ...]
if ~isa(M,'sym')
    convert_element = @(n) num2str(n);
else
    convert_element = @(n) char(n);
end

string = '';
for ii = 1:s1
    if ii==1; string=[string,'[']; end
    for jj=1:s2
        string = [string,convert_element(M(ii,jj)),','];
        if jj==s2; string = [string(1:end-1),';']; end
    end
    if ii==s1; string = [string(1:end-1),'];']; end
    if ii==1; fprintf(fid,'%s\n',string);
    else fprintf(fid,'\t%s\n',string);
    end
    string='';
end
fprintf(fid,'%s\n','  ');



end