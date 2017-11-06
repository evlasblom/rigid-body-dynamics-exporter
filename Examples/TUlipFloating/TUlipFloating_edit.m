function Model = TUlipFloating_edit(Model)

RigidBody = Model.rigidbody;

% Body Width
bw = 0.05;
Model.visual.bodywidth = bw;

% Virtual chian
% Force virtual chain to have zero size.
RigidBody(1).visual.bodysizemin = [0 0 0];
RigidBody(1).visual.bodysizemax = [0 0 0];
RigidBody(2).visual.bodysizemin = [0 0 0];
RigidBody(2).visual.bodysizemax = [0 0 0];
RigidBody(3).visual.bodysizemin = [0 0 0];
RigidBody(3).visual.bodysizemax = [0 0 0];
RigidBody(4).visual.bodysizemin = [0 0 0];
RigidBody(4).visual.bodysizemax = [0 0 0];
RigidBody(5).visual.bodysizemin = [0 0 0];
RigidBody(5).visual.bodysizemax = [0 0 0];

% Vertices and faces to make the trunk look fancy
w1 = 0.112;
w2 = 0.343;
h1 = 0.12;
h2 = 0.193;
h3 = 0.076;
th = 0.135;
[trunkVertices,trunkFaces] = ...
    make_trunk(w1,w2,h1,h2,h3,th,[0 -w1/2 -(h1+h2+h3)]);
RigidBody(6).visual.vertices = trunkVertices;
RigidBody(6).visual.faces = trunkFaces;

% Size of feet
% Final body lengths should be specified as they are not inside any of the
% homogeneous transformations.
fw = 0.15;
fl = 0.21;
off = 0.06;

RigidBody(12).visual.bodysizemin = -[off,fw/2,0.05];
RigidBody(12).visual.bodysizemax = [fl-off,fw/2,0];
RigidBody(18).visual.bodysizemin = -[off,fw/2,0.05];
RigidBody(18).visual.bodysizemax = [fl-off,fw/2,0];




% % Other
% imutrunk = r.IMUTrunk;
% dlowerleg = d.KneeYAnkleY;
% footheight = abs(r.AnkleXFoot(3));

% Output
Model.rigidbody = RigidBody;




% Vertices and faces to make a cube
    function varargout = make_trunk(w1,w2,h1,h2,h3,th,compos)
        % w1: width under
        % w2: width top
        % h1: height lower square
        % h2: height middle part
        % h2: height upper square
        % th: thickness
        wadd = (w2-w1)/2;
        hmid = h1+h2;
        htot = h1+h2+h3;
        coordinates = [0 0 0; 0 w1 0; th w1 0; th 0 0; ...
                       0 0 h1; 0 w1 h1; th w1 h1; th 0 h1 ; ...
                       0 -wadd hmid; th -wadd hmid; th w1+wadd hmid; 0 w1+wadd hmid ; ...
                       0 -wadd htot; th -wadd htot; th w1+wadd htot; 0 w1+wadd htot];
        coordinates(:,1) = coordinates(:,1) + compos(1);
        coordinates(:,2) = coordinates(:,2) + compos(2);
        coordinates(:,3) = coordinates(:,3) + compos(3);
        varargout{1} = coordinates;
        if nargout > 1 
            varargout{2} = [1 2 3 4; ...
                2 6 7 3; ... %
                4 3 7 8; ...
                1 5 8 4; ...
                1 2 6 5; ...
                6 12 11 7; ... %
                8 7 11 10; ...
                5 9 10 8; ...
                5 6 12 9; ...
                10 14 15 11; ... %
                12 11 15 16; ...
                9 13 16 12; ...
                9 10 14 13; ...
                13 14 15 16];
        end
    end
end