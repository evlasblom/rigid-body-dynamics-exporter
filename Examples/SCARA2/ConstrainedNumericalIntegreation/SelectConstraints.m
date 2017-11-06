function con_out = SelectConstraints(select,con_in)
% con_in:
% D     nc x 1 vector
% dD    nc x n matrix
% dconv nc x 1 vector
con_out = [];
cnt = 1;
for ii = 1:length(select);
    if select(ii) == 1;
        con_out(cnt,:) = con_in(ii,:);
        cnt =  cnt+1;
    end
end
end