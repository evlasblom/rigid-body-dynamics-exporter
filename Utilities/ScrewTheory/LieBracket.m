%LIEBRACKET Calculates the Lie bracket of two vector fields.
%   T = LIEBRACKET(T1,T2) calculates the twist resulting from the Lie 
%   bracket of two twists as T = T1*T2 - T2*T1.
%
%   Made by Erik Vlasblom
%   Last modified: 14-05-2013
function T = LieBracket(T1,T2)
if size(T1,1) == 4 && size(T1,2) && size(T2,1) && size(T2,2)
    T = T1*T2 - T2*T1;
elseif length(T1) == 6 && length(T2) == 6;
    Ta = Tilde(T1);
    Tb = Tilde(T2);
    T = UnTilde(Ta*Tb - Tb*Ta);
else
    error('The input should be a twist (4x4 or 6x1)')
end
end