function out = MTLB_write_atom(atom)
M = atom.exp;

% conversion method
if ~isa(M,'sym'); conv = @(n) num2str(n);
else conv = @(n) char(n);
end

% write expression
exp = ''; string = '';
for ii=1:size(M,1)
    if ii==1; string = [string,'[']; end

    for jj=1:size(M,2)
        string = [string,conv(M(ii,jj)),','];
        if jj==size(M,2); string = [string(1:end-1),';']; end
    end
    
    if ii==size(M,1); 
        string = [string(1:end-1),']'];
        newline = '';
    else
        newline = '\n';
    end
        
    if ii==1; exp = sprintf(['%s%s',newline],exp,string);
    else exp = sprintf(['%s\t%s',newline],exp,string);
    end
    
    string='';
end

out.name = atom.name;
out.num = atom.num;
out.caller = atom.caller;
out.exp = exp;

end