%ROTATIONMATRIX Calculates the desired rotation matrix in 3D.
%   R = ROTATIONMATRIX(A,o) gives the rotation of angles A in the order o.
%   The angles, in radians, can either be separated by comma's or put in a 
%   vector. The order should be a string. Any rotation order can be used
%   for 1, 2 or 3 axes. The rotation matrix rotates from a local frame x' 
%   to a global frame x: x = R x'. Use the inverse or transpose for the
%   ohter way around.
%  
%   Example:
%   A = [45,90,90]*pi/180;
%   R = matrix2la(A,'zxz')
%
%   R =
%          0   -0.7071    0.7071
%          0   -0.7071   -0.7071
%     1.0000         0         0
%
%   Which rotates 45 degrees about z, 90 about x and 90 about z.
%   
%   Made by Pepijn Cox, Erik Vlasblom
%   Last modified: 17-04-2013
function R = RotationMatrix(varargin)
%------------------------------------
% Import inputs in the right vars
ninputs = nargin;
useInverse = false;
if nargin > 0
    if size(varargin{end},1)==1 && size(varargin{end},2)==1 && isnumeric(varargin{end})
        ninputs = nargin - 1;
        if varargin{end} == -1
            useInverse = true;
        end
    end
end

switch ninputs
    % no a
    case 0
        e = 0;
        typeRotation = 'n';
    % a is a vector
    case 2
        if not(length(varargin{1}) == length(varargin{2}))
            error('The lengths of the input angles and rotation order do not match.');
        end
        if(length(varargin{1}) == 3) || (length(varargin{1}) == 2) || ...
                (length(varargin{1}) == 1)
            e = varargin{1};
            typeRotation = lower(varargin{2});
        else
            error('Wrong input was given.');
        end
    otherwise
        error('Wrong input was given.');
end


%------------------------------------
% Calculating individual rotations
tmpsin = sin(e);
tmpcos = cos(e);
    function rot = get_r(axis,a)
        switch axis
            case 'x'
                rot = [1              0               0
                    0              tmpcos(a)      -tmpsin(a)
                    0              tmpsin(a)      tmpcos(a)];
            case 'y'
                rot = [tmpcos(a)     0               tmpsin(a)
                    0              1               0
                    -tmpsin(a)    0               tmpcos(a)];
            case 'z'
                rot = [tmpcos(a)     -tmpsin(a)     0
                    tmpsin(a)     tmpcos(a)      0
                    0              0               1];
        end
    end


%------------------------------------
% Calculating total rotation
r = eye(3);

for i=1:length(e)
    switch typeRotation(i)
        case 'x'
            x = get_r('x',i);
            r = r*x;
        case 'y'
            y = get_r('y',i);
            r = r*y;
        case 'z'
            z = get_r('z',i);
            r = r*z;
        case 'n'
            r = eye(3);
        otherwise
            error('Error: Invalid input Euler angle order type (conversion string)');
    end
end

R = r;

%------------------------------------
% Give output

if ~useInverse
    R = R;
else
    R = R.';
end

end