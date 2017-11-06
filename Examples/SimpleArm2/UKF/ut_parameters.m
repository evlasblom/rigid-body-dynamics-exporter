function ut = ut_parameters(varargin)

if nargin > 0
    set = varargin{1};
    L = varargin{2};
else
    set = 'default';
end

ut.set = set;

switch set
    case 'default'
        ut.k = 0;
        ut.a = 1e-3;
        ut.b = 2;
        ut.lam = ut.a^2*(L+ut.k)-L;
        ut.weight = 1/(L+ut.lam);
        
        ut.scale = L+ut.lam;
        ut.Wm = [ut.weight*ut.lam 0.5*ut.weight*ones(1,2*L)];
        ut.Wc = [ut.Wm(:,1)+(1-ut.a^2+ut.b) ut.Wm(:,2:end)];
    case 'general_symmetric'
        ut.scale = L;
        ut.Wm = [0 1/(2*L)*ones(1,2*L)];
        ut.Wc = ut.Wm;
    case 'simplex'
        ...
    otherwise
        warning('No existing simplex set chosen')
end

if ~isfield(ut,'Wm') || ~isfield(ut,'Wc')
    error('No weighting was set for the unscented transformation')
elseif ~isfield(ut,'scale')
    error('No scaling was set for the cholesky factor')
end