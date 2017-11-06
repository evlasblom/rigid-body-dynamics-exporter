%COMCalculator Calculates the Center of Mass.
%   r = COMCALCULATOR(model,q,parameters) computes the Center of Mass of 
%   the model using the provided joint angles q. The parameter structure
%   can be used if some parameters require changes. If not provided, the
%   default values from the model are used.
%
%   Inputs:
%   model       model structure as created by CreateModel
%   q           joint angles
%   parameters  parameter structure with modified parameters (optional)

%   Made by Erik Vlasblom
%   Last modified: 21-06-2015
function r = COMCalculator(model,q,varargin)

% modify parameters
if nargin > 2
    par = varargin{1};
    if isa(par,'struct')
        model = umfp(model,par);
    else
        error('Wrong input')
    end
end

% pre-allocate
Hglobal = eye(4);
m = 0;
r = zeros(3,1);

% compute com position
B = model.rigidbody;
for ii = 1:model.dof
    Hlocal = CreateTransform(B(ii).joint.type,B(ii).joint.axis,B(ii).joint.offset,q(ii));
    Hglobal = Hglobal*Hlocal;
    Hm = HomogeneousTransform(B(ii).inertial.origin);
    r_ii = h2p(Hglobal*Hm);
    
    m = m + B(ii).inertial.mass;
    r = r + B(ii).inertial.mass*r_ii;
end
r = r./m;

end