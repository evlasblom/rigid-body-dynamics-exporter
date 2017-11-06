function atom = MTLB_atom(name,num,value)
% ATOM = {name,numbers,caller,expression}
atom.name = name;
atom.num = num;
atom.caller = MTLB_name(atom.name,atom.num);
atom.exp = value; % numeric / symbolic
atom.type = 'atom';
atom.zero = false;

end