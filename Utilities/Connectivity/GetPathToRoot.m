%GETPATHTOROOT gets the joints from body i to the root.
%   kappa = GETPATHTOROOT(lambda,i) uses the connectivity to obtain the
%   joints that lie in between body i and the root, including body i,
%   exluding the root.
%
%   Made by Erik Vlasblom
%   Last modified: 12-09-2014
function kappa = GetPathToRoot(lambda,i)
if i > 1
    kappa = [i, GetPathToRoot(lambda,lambda{i})];
else
    kappa = 1;
end
end