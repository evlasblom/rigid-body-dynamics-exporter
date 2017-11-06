%INVERSEHOMOGENEOUS inverts a homogeneous transformation.
%   Hinv = INVERSEHOMOGENEOUS(H) calculates the inverse of homogeneous transformation H.
%     
%   Erik Vlasblom
%   Last modified: 17-01-2014
function Hinv = InverseHomogeneous(H)
%------------------------------------
R = H(1:3,1:3);

if (size(H,1) == 4) && (size(H,2) == 4)
    p = H(1:3,4);
    Hinv = [R.' -R.'*p; zeros(1,3) 1];
elseif (size(H,1) == 3) && (size(H,2) == 3)
    Hinv = R.';
else
    error('Input should be a 4x4 transformation or 3x3 rotation matrix')
end

end