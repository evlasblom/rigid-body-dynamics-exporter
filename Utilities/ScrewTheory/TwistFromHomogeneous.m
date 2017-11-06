%TWISTFROMHOMOGENEOUS Obtains the twist from a homogeneous transformation.
%   T = TWISTFROMHOMOGENEOUS(H,theta) computes the relative twist from the
%   homogeneous transformation H such that H(theta) = e^(T*theta) * H(0).
%   For a transformation H a body frame b to a spatial frame a, the
%   calculated twist is the local twist that describes the velocities of
%   the body relative to the spatial frame, expressed in the spatial frame.
%
%   T = TWISTFROMHOMOGENEOUS(H,theta,1) computes the relative twist as
%   above, but expressed in the body frame.
%
%   Made by Erik Vlasblom
%   Last modified: 07-05-2013
function T = TwistFromHomogeneous(H,theta,varargin)
if isa(H,'sym')
    dH = diff(H,theta);
    T = simplify(dH*InverseHomogeneous(H)); % spatial
    if nargin>2
        T = simplify(InverseHomogeneous(H)*dH); % body
    end
else
    error('The input is not symbolic')
end
end