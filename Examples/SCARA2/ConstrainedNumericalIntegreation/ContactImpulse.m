function xo = ContactImpulse(model,FDYN,FCON,e,xn,imp,n,lam)
q = xn(1:n);
dq = xn(n+1:end);
[M,~,~] = FDYN(xn);
[~,dD,~] = FCON(xn);
%ctot = sum(imp);
xo = xn;
vminus = dq;
dD = SelectConstraints(imp,dD);
Mt = [M, -dD.'; dD, zeros(sum(imp))];
mom = [M*vminus; -e*dD*vminus];
%Mi = Mt'/(Mt*Mt' + lam^2*eye(n+ctot));
%a = Mi*mom;
a = Mt\mom;
xo(n+1:end) = a(1:n);
end