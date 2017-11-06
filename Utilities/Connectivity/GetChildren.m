%GETCHILDREN gets the joints from body i to the root.
%   mu = GETCHILDREN(lambda,i) uses the connectivity to obtain the
%   children bodies of body i.
%
%   Made by Erik Vlasblom
%   Last modified: 12-09-2014
function mu = GetChildren(lambda,i)
c = 0;
mu = cell(0);
for ii = 1:length(lambda)
    if lambda{ii} == i
        c = c+1;
        mu{c} = ii;
    end
end
end