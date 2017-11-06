%H_TRANSFORM Calculates a transformation matrix of a kinematic chain.
%   H = H_TRANSFORM(model,q,from,to) provides the 
%   homogeneous transformation matrix to express a point given in a body's 
%   local coordinate system in the coordinate system of another body
%   according to the kinematic structure of the model
%
%   Inputs:
%   model       model structure as created by CreateModel
%   q           joint angles
%   from        frame id of the current coordinate system
%   to          frame id in which the transformed point should be expressed
%
%   Made by Erik Vlasblom
%   Last modified: 27-05-2015
function H = H_Transform(model,q,from,to)

%------------------------------------
% Find shortest path from local to target using Dijkstra's algorithm

% number of bodies includes the base
l = length(model.connectivity.lambda)+1;

% append a -1 to the path to mark the end of the path
path = -ones(1,1,l+1); 

% caculate path
if to == 0 && from == 0
    path(1) = 0;
elseif to == 0
    toroot = [model.connectivity.kappa{from}, 0];
    path(1:length(toroot)) = fliplr(toroot);
elseif from == 0
    toroot = [model.connectivity.kappa{to}, 0];
    path(1:length(toroot)) = toroot;
else
    % add 1, do dijkstra, subtract 1 (because 0 is not a valid index)
    path(1:l) = Dijkstra(model.connectivity.G,[from+1,to+1])-1;
end


%------------------------------------
% Loop through the path to find the total transformation
H = eye(4);
B = model.rigidbody;

ii = 1;
from = path(ii);
to = path(ii+1);
while to ~= -1
    if to == 0;
        to_parent = -1;
    else
        to_parent = model.connectivity.lambda{to};
    end
    
    if to_parent == from % if next body is parent
        Hlocal = CreateTransform(B(to).joint.type,B(to).joint.axis,B(to).joint.offset,q(to));
    else % if next body is child
        Hlocal = CreateTransform(B(from).joint.type,B(from).joint.axis,B(from).joint.offset,q(from));
        Hlocal = InverseHomogeneous(Hlocal);
    end
    
    H = H*Hlocal;
    
    ii = ii+1;
    from = path(ii);
    to = path(ii+1);
end


end
