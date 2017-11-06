%TRANSFORM Translates and rotates vectors with x-y-z coordinates.
%   Dt = TRANSFORM(D,H) uses the homogeneous transformation H to transform 
%   the 3-by-n or n-by-3 matrix D. Note that H is 4-by-4.
%
%   Dt = TRANSFORM(D,a,o,p) creates a homogeneous transformation matrix
%   using the HomogeneousTransform(a,o,p) and applies this to all coordinates that
%   reside in the 3-by-n or n-by-3 D matrix. The angles a and the
%   translations p can either be vectors or separate arguments.
%
%   Made by Erik Vlasblom
%   Last modified: 02-05-2013
function Dnew = Transform(D,varargin)

%------------------------------------
% Import inputs in the right vars
switch nargin
    case 2
        if size(varargin{1},1) == 4 && size(varargin{1},2) == 4
            H = varargin{1};
        else
            H = HomogeneousTransform(varargin{1});
        end
    case 3
        H = HomogeneousTransform(varargin{1},varargin{2});
    case 4
        H = HomogeneousTransform(varargin{1},varargin{2},varargin{3});
    case 5
        H = HomogeneousTransform(varargin{1},varargin{2},varargin{3}, ...
            varargin{4});
    otherwise
        error('Wrong input was given.');
end

rows = size(D,1);
if rows ~= 3
    D = D.';
end

%------------------------------------
% Do transformation and give the right output
Dtrans = H*[D; ones(1,size(D,2))];

if rows ~= 3
    Dnew = Dtrans(1:3,:).';
else
    Dnew = Dtrans(1:3,:);
end

end