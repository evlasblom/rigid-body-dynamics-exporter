function out = MTLB_coriolis_element(name,num,dM,dq,n)
% The implementation of the coriolis elements is somewhat dirty. It is not
% following the same logic as all other functions.

out.name = name;
out.num = num;
out.caller = [name,'(',regexprep(num2str(num),'\s+',','),')'];

i = num(1);
j = num(2);
dMi = dM{i};
dMj = dM{j};

zero = true;
string = '';
for k = 1:n
    dMk = dM{k}; dqk = dq(k);
    
    if dMk.zero && dMj.zero && dMi.zero
        string = '0 + ';
    else     
        zero = false;
        
        string = [string,'0.5*('];
        if ~dMk.zero
            string = [string,'+',dMk.caller,'(',num2str(i),',',num2str(j),')'];
        end
        if ~dMj.zero
            string = [string,'+',dMj.caller,'(',num2str(i),',',num2str(k),')'];
        end
        if ~dMi.zero
            string = [string,'-',dMi.caller,'(',num2str(k),',',num2str(j),')'];
        end
        string = [string,')*',char(dqk),' + '];
    end
end

out.exp = string(1:end-2);

out.zero = zero;

end