function [l,act,stab,loc,imp,returnflag] = ConstraintDetection(model,FDYN,FCON,FCONDET,beta,xn,un,act,stab,loc,imp,n,tol,lam)
returnflag = '';
[Dnum,~,~] = FCON(xn);

if sum(act) > 0
    [~,l] = ForwardDynamics(...
        model,...
        FDYN,...
        FCON,...
        beta,xn,un,act,n,lam);
else
    l = zeros(1,length(act));
end

[act,stab,loc,imp] = FCONDET(model,xn,Dnum,l,act,stab,loc,imp,tol);

end