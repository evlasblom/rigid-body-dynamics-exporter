%INITIALCONFIGURATION Returns te initial configuration of a transformation.
%   H0 = INITIALCONFIGURATION(H, theta) returns the initial configuration of a
%   transform H such that H = TwistExponent((TwistFromHomogeneous(H))*H0.
%
%   Made by Erik Vlasblom
%   Last modified: 05-02-2013
function H0 = InitialConfiguration(H,theta)
H0 = subs(H,theta,zeros(size(theta)));
end