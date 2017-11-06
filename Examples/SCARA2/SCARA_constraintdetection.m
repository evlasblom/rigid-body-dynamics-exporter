function [act,stab,loc,imp,returnflag] = SCARA_constraintdetection(model,xn,Dnum,l,act,stab,loc,imp,tol)

% Constraint 1, maximum angle theta2
cc = 1;
if act(cc) == 0;
    if Dnum(cc) <= -tol;
        act(cc) = 1; stab = act; imp = act;
        loc(cc) = 0;
    end
elseif act(cc) == 1;
    if l(cc) <= 0
        loc(cc) = 0;
        act(cc) = 0;
    end
end

% Constraint 2, friction, body 4
% Constraint 3, friction, body 4
% Constraint 4, touchdown, body 4
mu = 0.6; 
cc = 4;
cca = [2 3 4];
if act(cc) == 0;
    if Dnum(cc)-0.05 <= -tol;
        act(cca) = 1; stab = act; imp = act; 
        H4_0 = H_Transform(model,xn,4,0);
        pcontact = Transform([0,0,-0.05],H4_0);
        loc(cca) = [pcontact(1) pcontact(2) 0]; 
    end
elseif act(cc) == 1;
    if l(cc) <= 0
        loc(cca) = 0;
        act(cca) = 0; 
        returnflag = 'release';
    elseif sqrt(l(cca(1))^2+l(cca(2))^2) >= mu*l(cc)
        loc(cca) = 0;
        act(cca) = 0;
        returnflag = 'slipping';
    end
end

end
