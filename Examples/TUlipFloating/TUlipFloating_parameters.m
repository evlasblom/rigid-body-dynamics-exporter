%TULIP_PARAMETERS Defines the required parameters values (pv).
%   [d,r,m,i] = TULIP_PARAMETERS(pn) sets the parameters of TUlip that are
%   required for the derivation of the equations of motion. The true
%   parameters are modified such that the FFH model matches TUlip.
%
%   d: offset between joints
%   r: center of mass offset from body frame
%   m: mass
%   i: inertia
%
%   Main source for parameters:
%   https://robocup.3xo.eu/trac/ROBOCUP-3TU/ROBOCUP2008/browser
%
%   Made by Erik Vlasblom
%   Last modified: 15-09-2014

function [dout,rout,mout,iout] = TUlipFloating_parameters()

% found an error! 
% upperleg origin should be at the top in the center of the leg
% but it is at the top in the inside of the leg
origin_correction = 0.02;

bodynames = {'Trunk'; 'HipXMotor'; 'HipYAxle'; 'UpperLeg'; 'LowerLeg'; 'AnkleYAxle'; 'Foot'};

% BODY LENGTHS [m]
% The values are found for the left leg.
% obtained from /Tulip/trunk/software/TUlipMC/src/platform/kinematics/TUlipKinematicModelFactory.hpp
% on 24-07-2013
d.IMUTrunk = [0.002 0 0.242]; %
d.TrunkHipZ = [0.085 0.07594 -0.392]; %
d.HipZHipX = [0 0 -0.0705]; 
d.HipXHipY = [0 0.03775+origin_correction -0.04]; 
d.HipYKneeY = [0 0 -0.275];
d.KneeYAnkleY = [0 0 -0.32];
d.AnkleYAnkleX = [0 0 0];
d.AnkleXFoot = [0 0 -0.05];


% CENTER OF MASS POSITIONS [m]
% Defined in body frame that has the same origin as the parent joint
% The values are found for the right leg.
% obtained from /Tulip/trunk/software/TUlipMC/src/platform/kinematics/TUlipKinematicModelFactory.hpp
% on 24-07-2013
r.TrunkMeas = [-0.0235 0.09656 -0.427];
r.TrunkCOMRel = [0.099 -0.0975 0.225];
r.HipXMotorMeas = [0.045 0 0]; 
r.HipXMotorCOMRel = [-0.054 0 0.007]; % --> [0 0 0]
r.HipYAxleMeas = [0 0 0];
r.HipYAxleCOMRel = [0 0 0];
r.UpperLegMeas = [-0.055 0 -0.015];
r.UpperLegCOMRel = [0.05369 -0.0302+origin_correction -0.08228];
r.LowerLegMeas = [0 -0.04+origin_correction 0];
r.LowerLegCOMRel = [0.0375 0.02059 -0.155];
r.AnkleYAxleMeas = [0 0 0];
r.AnkleYAxleCOMRel = [0 0 0];
r.FootMeas = [-0.06 0.063 -0.001];
r.FootCOMRel = [0.08714 -0.06848 -0.00454];
for i = 1:length(bodynames);
    r.(bodynames{i}) = r.([bodynames{i},'Meas']) + r.([bodynames{i},'COMRel']);
end
r.Trunk(1) = 0.06;
r.UpperLeg(1) = 0.005;
r.LowerLeg(3) = -0.15;
r.Foot = [0.035 -0.01 -0.01];


% MASSES [kg]
% obtained from /Tulip/trunk/software/TUlipMC/src/platform/kinematics/TUlipConstants.h
% on 24-07-2013
massTrunk = 8.1540;
massHipXMotor = 0.61;
massHipYAxle = 0.075;
massUpperLeg = 2.14;
massLowerLeg = 1.023;
massAnkleYAxle = 0.075;
massFoot = 0.37;

% INERTIAS [kg*m^2]
% obtained from /Tulip/branches/software/branch_TUe/TUlipMC/src/statemachine/behaviours/TulipParameters.h
% on 24-07-2013
inertiaTrunk = [0.27217, 0.20518, 0.07047];
inertiaHipXMotor = [0.00002, 0.00095, 0.00056]; % --> 0 0 0
inertiaHipYAxle = [0.00002, 0.00095, 0.00056]; % --> 0 0 0
inertiaUpperLeg = [0.01467, 0.01134, 0.00148]; 
inertiaLowerLeg = [0.00396, 0.00368, 0.00009]; 
inertiaAnkleYAxle = [0.00002, 0.00095, 0.00056]; % --> 0 0 0
inertiaFoot = [0.00040, 0.00157, 0.00205];


% OUTPUT
dout = [zeros(6,3); % floating body
        d.TrunkHipZ; % left leg
        d.HipZHipX; 
        d.HipXHipY; 
        d.HipYKneeY;
        d.KneeYAnkleY;
        d.AnkleYAnkleX;
        mirror_y(d.TrunkHipZ); % right leg
        mirror_y(d.HipZHipX);
        mirror_y(d.HipXHipY); 
        mirror_y(d.HipYKneeY)
        mirror_y(d.KneeYAnkleY);
        mirror_y(d.AnkleYAnkleX)];
    
rout = [zeros(5,3) % floating body
        r.Trunk; 
        mirror_y(r.HipXMotor); % left leg 
        mirror_y(r.HipYAxle);
        mirror_y(r.UpperLeg);
        mirror_y(r.LowerLeg);
        mirror_y(r.AnkleYAxle);
        mirror_y(r.Foot);
        r.HipXMotor; % right leg 
        r.HipYAxle;
        r.UpperLeg;
        r.LowerLeg;
        r.AnkleYAxle;
        r.Foot];
    
mout = [zeros(5,1); % floating body
        massTrunk; 
        massHipXMotor; % left leg
        massHipYAxle;
        massUpperLeg;
        massLowerLeg; 
        massAnkleYAxle;
        massFoot;
        massHipXMotor; % right leg
        massHipYAxle;
        massUpperLeg;
        massLowerLeg; 
        massAnkleYAxle;
        massFoot];
    
iout = [zeros(5,3); % floating body
        inertiaTrunk; 
        inertiaHipXMotor; % left leg
        inertiaHipYAxle; 
        inertiaUpperLeg;
        inertiaLowerLeg;
        inertiaAnkleYAxle;
        inertiaFoot;
        inertiaHipXMotor; % right leg
        inertiaHipYAxle; 
        inertiaUpperLeg;
        inertiaLowerLeg;
        inertiaAnkleYAxle;
        inertiaFoot];
    
% % OTHER OUTPUT (use later)
% % parameters for the trunk (for animation)
% pv.o.imutrunk = r.IMUTrunk;
% pv.o.dlowerleg = d.KneeYAnkleY;
% pv.o.w1 = 0.112;
% pv.o.w2 = 0.343;
% pv.o.h1 = 0.12;
% pv.o.h2 = 0.193;
% pv.o.h3 = 0.076;
% pv.o.th = 0.135;
% % foot height
% pv.o.footheight = abs(r.AnkleXFoot(3));

% --------------------------------------------------
% Small function to mirror values over the xz-plane.
    function yout = mirror_y(yin)
        yout = yin;
        yout(2) = -yout(2);
    end

end