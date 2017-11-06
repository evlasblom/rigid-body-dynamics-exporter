function MTLB_export_all(options,list,x,parameters)

% Initialize
WR_ATOM = @MTLB_write_atom;
WR_STATE = @MTLB_write_statestring;
WR_PAR = @MTLB_write_parameterstring;

str = '';

datanames = '';
for ii = 1:length(options.exportdata);
    datanames = [datanames,',',options.exportdata{ii}];
end
datanames = datanames(2:end);
f_name = [options.name,'_',regexprep(datanames,',','_')];

fprintf(['  - function: ',f_name,'\n']);

% Header
f_output = ['[',datanames,']'];
f_input = 'q';
if options.par_edit; f_input = [f_input,',p']; end
headstring = sprintf('%s\n\n\n',['function ',f_output,' = ',f_name,'(',f_input,')']);
statestring = WR_STATE(x);
parstring = '';
if options.par_edit; parstring = WR_PAR(parameters); end
str = sprintf('%s%s%s%s\n',str,headstring,statestring,parstring);

% Content
for ii = 1:length(list)
    if isfield(list{ii},'type')
        if strcmp(list{ii}.type,'atom')
            list{ii} = WR_ATOM(list{ii});
        end
    end
    lines = sprintf('%s%s%s%s\n',list{ii}.caller,' = ',list{ii}.exp,';');
    str = sprintf('%s\n%s',str,lines);
end
str = sprintf('%s\n\n',str);

% Sub functions (matlab)
str = sprintf('%s\n\n%s\n\n',str,'% ---------- sub functions');
func = {'f_adjoint.m','f_inverseadjoint.m','f_ad.m','f_skew.m','f_rbar.m','f_ombar.m'};
for ii = 1:length(func)
    fid = fopen(func{ii});
    tline = fgetl(fid);
    first = true;
    while ischar(tline)
        if first
            str = sprintf('%s\t%s\n',str,tline);
            first = false;
        elseif strcmp(tline,'end')
            str = sprintf('%s\t%s\n\n',str,tline);
        else
            str = sprintf('%s\t\t%s\n',str,tline);
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

% Footer
str = sprintf('%s\n%s',str,'end');

% Print to file
fid = fopen([options.exportloc,f_name,'.m'],'w');
fprintf(fid,'%s',str);
fclose(fid);

end