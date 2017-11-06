%CREATETRANSFORM Create the specified homogeneous transform.
%   H = CREATETRANSFORM(t,a,d,q) creates a homogeneous transformation 
%   matrix for a body with a specific joint type. It uses properties from 
%   the model such as the joint type (t), the axis of rotation or 
%   translation (a) and the constant offset (d). The degree of freedom 
%   should be specified by a symbolic variable (q).
function H = CreateTransform(type,axis,d,q)

if type == jointtype.Revolute
    H = HomogeneousTransform(q,axis,d);
elseif type == jointtype.Prismatic
    ax = ('xyz' == axis);
    H = HomogeneousTransform(ax*q+d);
else
    error('No such jointType');
end

end