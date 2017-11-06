% Converts model structure for RBDE (1.x.x) to the model structure that is
% required for Spatial v2.
%
% Not yet checked for errors!
% Namely the use of spatial definitions... (e.g. xlt function)
function spatial_model = rbde2spatial(rbde_model)

mdl = rbde_model;

% ==========
% Notes on compatibility between Screw Theory and Spatial Algebra:
%
% a^X_b = Ad2X(InverseAdjoint(H^b_a))
% with a^X_b(r) and H^b_a(r) using the same r.
% with Ad2X = @(Ad) [Ad(1:3,1:3) Ad(4:6,1:3); Ad(1:3,4:6) Ad(4:6,4:6)]
% and I2Is = @(I) [I(4:6,4:6), I(4:6,1:3).'; I(1:3,4:6).', I(1:3,1:3)] 
%
% ==========
Ad2X = @(Ad) [Ad(1:3,1:3) Ad(4:6,1:3); Ad(1:3,4:6) Ad(4:6,4:6)];

% --- OK Dof, lambda
new.NB = mdl.dof;
new.parent = cell2mat(mdl.connectivity.lambda);
mu = mdl.connectivity.mu;

% --- OK Joint type
new.jtype = cell(1,new.NB);
for ii =  1:new.NB
    if mdl.rigidbody(ii).joint.type == jointtype.Prismatic;
        j = 'P';
    elseif mdl.rigidbody(ii).joint.type == jointtype.Revolute;
        j = 'R';
    else
        error('No such joint')
    end
    new.jtype{ii} = [j,mdl.rigidbody(ii).joint.axis];
end

% --- OK Inverse adjoint
new.Xtree = cell(1,new.NB);
for ii = 1:new.NB
    new.Xtree{ii} = xlt(mdl.rigidbody(ii).joint.offset);
    %new.Xtree{ii} = Ad2X(InverseAdjoint(HomogeneousTransform(mdl.rigidbody(ii).joint.offset)));
end

% --- OK Inertia at body frame
new.I = cell(1,new.NB);
for ii = 1:new.NB
    m = mdl.rigidbody(ii).inertial.mass;
    I = mdl.rigidbody(ii).inertial.i;
    r = mdl.rigidbody(ii).inertial.origin; 
    new.I{ii} = xlt(r).'*diag([I m m m])*xlt(r);
    %tmpAd = @(r) Ad2X(InverseAdjoint(HomogeneousTransform(r)));
    %new.I{ii} = tmpAd(r).'*diag([I m m m])*tmpAd(r);
end

% --- Gravity
new.gravity = mdl.gravity;

% --- Joint constraints
% - gamma_q not used

% --- OK Visual
% Obtain body lengths from parameters (copied from Animation3D.m)
B = mdl.rigidbody;
for k = 1:new.NB
    % Obtain body length from transformation to the body's children
    bw = mdl.visual.bodywidth;
    maxd = [bw, bw, bw]/2;
    mind = -[bw, bw, bw]/2;
    for ii = 1:length(mu{k});
        maxd = max(maxd,B(mu{k}{ii}).joint.offset);
        mind = min(mind,B(mu{k}{ii}).joint.offset);
    end
    % Use body properties that are set, else, find body size
    if isempty(B(k).visual.bodysizemax);
        B(k).visual.bodysizemax = maxd;
    end
    if isempty(B(k).visual.bodysizemin)
        B(k).visual.bodysizemin = mind;
    end
end
mdl.rigidbody = B;
% Set for spatial
new.appearance.body = cell(1,new.NB);
for ii = 1:new.NB
    vertices = [mdl.rigidbody(ii).visual.bodysizemin;
                mdl.rigidbody(ii).visual.bodysizemax];
    if all(all(vertices == 0)) == true 
        new.appearance.body{ii} = {};
    else
        r = mdl.rigidbody(ii).inertial.origin;
        comsphere = { 'facets', 80, ...
            'colour', [0.1 0.1 0.1; 0.8 0.8 0.8], ...
            'sphere', r, bw };
        new.appearance.body{ii} = {'box', vertices,...
            comsphere};
    end
end

% --- Camera view
% - camera not used

% --- Ground contact
% - gc not used

spatial_model = new;

end

