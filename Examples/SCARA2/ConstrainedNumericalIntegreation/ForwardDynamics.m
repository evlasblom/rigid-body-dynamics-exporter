function [dx,l] = ForwardDynamics(model,FDYN,FCON,beta,x,u,act,n,lam)

q = x(1:n);
dq = x(n+1:end);
[M,C,G] = FDYN(x);
[~,dD,dconv] = FCON(x);
ctot = sum(act);
if ctot > 0
    dD = SelectConstraints(act,dD);
    dconv = SelectConstraints(act,dconv);
    Mt = [M, -dD.'; dD, zeros(sum(act))];
    ft = [-C*dq-G-beta*dq+u; -dconv];
else
    Mt = M;
    ft = -C*dq-G-beta*dq+u;
end
% Mi = Mt'/(Mt*Mt' + lam^2*eye(n+ctot));
% ddq = Mi*ft;
ddq = Mt\ft;
dx = [dq; ddq(1:n)];
if nargout > 1
    l = [];
    if ctot > 0
        l = zeros(1,length(act));
        l(1,act==1) = ddq(n+1:end);
    end
end

end