%UMFP update model from parameters
%   Model = UMFP(Model,par) uses the values of the parameter structure to
%   update the Model structure.
%
%   Made by Erik Vlasblom
%   Last modified: 10-10-2014
function Model = umfp(Model,parameters)

pnames = fieldnames(parameters);
for ii = 1:length(pnames)
    pname = pnames{ii};
    [numloc1,numloc2] = regexp(pname,'\d*'); % location of number
    
    par = pname(1:numloc1-1); % type of parameter
    num = str2double(pname(numloc1:numloc2)); % body number
    if length(pname) == numloc2; el = 1; % element (x y or z)
    else el = regexpi('xyz',pname(end));
    end
    
    switch par
        case 'd'
            Model.rigidbody(num).joint.offset(el) = parameters.(pname);
        case 'r'
            Model.rigidbody(num).inertial.origin(el) = parameters.(pname);
        case 'm'
            Model.rigidbody(num).inertial.mass(el) = parameters.(pname);
        case 'i'
            Model.rigidbody(num).inertial.i(el) = parameters.(pname);
        otherwise
                error('No such parameter')
    end
end
end