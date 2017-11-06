function xo = GaussNewton(FCON,xn,stab,loc,tol,lam)
%tol = 1e-12;
maxit = 5;
xh = xn(1:end/2);
%ctot = sum(stab);
[Dc,dD,~] = FCON(xn);
Dc = SelectConstraints(stab,Dc);
Dc = Dc - SelectConstraints(stab,loc);
dD = SelectConstraints(stab,dD);
j = 1;
while (max(abs(Dc)) > tol) && (j < maxit)
    % calculate: dx = -dD'*inv((dD*dD'))*D;
    DD = dD*dD';
    %DDi = DD'/(DD*DD' + lam^2*eye(ctot));
    %DDiD = DDi*Dc;
    DDiD = DD\Dc;
    dx = -dD'*DDiD;
    xh = xh + dx;
    xn(1:end/2) = xh;
    [Dc,dD,~] = FCON(xn);
    Dc = SelectConstraints(stab,Dc);
    Dc = Dc - SelectConstraints(stab,loc);
    dD = SelectConstraints(stab,dD);
    if isempty(Dc)
        Dc = tol/10;
        break
    end
    j = j + 1;
end
xo = xn;
end