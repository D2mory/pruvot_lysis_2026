% script virusremoval_data.m
%
% Generates data for the study of virus removal as a means to quantify
% the role of marine viruse in carbon cycling;
% see Figs. 4, 5, S2, S3 and S4 in Pruvôt et al. (2026)
%
% Data is saved in virusremoval_data.mat

% Fix parameter ranges
pr=struct('qP',[3e-7 3e-6], ...
          'qB',[2e-9 2e-8], ...
          'qG',[2e-4 2e-3], ...
          'qV',[2e-12 6e-11], ...
          'muP_max',[.2 2], ...
          'KP',[1e7 1e9], ...
          'lP_inorg',[.001 .1], ...
          'lP_org',[.005 .1], ...
          'psiPG',[1e-6 2e-5], ...
          'muB_max',[.5 2], ...
          'KDB',[30 300], ...
          'pDB_assim',[.05 .2], ...
          'lB_inorg',[.001 .1], ...
          'lB_org',[.005 .1], ...
          'psiBG',[1e-6 3e-5], ...
          'lG_inorg',[1e-3 1e-2], ...
          'lG_cons',[1e-7 1e-5],...
          'pG_assim',[.3 .5], ...
          'pG_org',[.2 .4], ...
          'phiP',[1e-9 1e-8], ...
          'betaP',[50 500], ...
          'lVP',[.05 5], ...
          'phiB',[1e-11 1e-9], ...
          'betaB',[10 100], ...
          'lVB',[.05 5]);

% Fix variable ranges
vr=[2e5 2e7;   % P
    2e7 2e9;   % B
    2e3 2e5;   % G
    2e6 2e8;   % VP
    2e8 2e10;  % VB
    30 100];   % D

%% Generate 10000 parameter sets with plausible equilibrium state

nSim = 1e4;
pars(nSim) = pr;            % structure to store selected parameter sets

fprintf('Starting simulations of full ecosystem model\n');

c1 = 1;
c2 = 1;
while c2 <= nSim
    pn = virusremoval_model_smp(pr);
    xeq = virusremoval_model_eqs(pn);
    if virusremoval_model_chk(xeq,pn,vr)
        if mod(c2,100) == 0
            fprintf('. found parameter set %d (out of 10000) at simulation %d\n', c2, c1)
        end
        pars(c2) = pn;
        c2 = c2+1;
    end
    c1 = c1+1;
end

%% Remove viruses and compute new equilibrium state (if any)

solref=struct('with',[],'withC',[],'wout',[],'woutC',[]);
% structure to store simulation results
% 'with' : solution with viruses, abundances
% 'withC' : solution with viruses, carbon content
% 'wout' : solution without viruses, abundances
% 'woutC' : solution without viruses, carbon content
sols(nSim)=solref;

fprintf('Starting simulations of virus removal\n');

nrms=zeros(1,nSim);
for c1=1:nSim
    if mod(c1,100) == 0
        fprintf('. found new steady state for parameter set %d (out of 10000)\n', c1)
    end

    %%% With viruses %%%
    pn = pars(c1);
    xeq = virusremoval_model_eqs(pn);
    convC = [pn.qP pn.qB pn.qG pn.qV pn.qV 1]';
    xeqC = convC.*xeq;
    sols(c1).with=xeq;
    sols(c1).withC=xeqC;

    %%% Without viruses %%%
    x0C = xeqC;
    x0C(4:5) = 0;              % virus removal in initial conditions
    opts = odeset('abstol',1e-8,'reltol',1e-8,'nonnegative',1:6); 
    [TC,XC] = ode45(@(t,xc) virusremoval_model_rhsC(xc,pn),[0 1e4],x0C,opts);
    xwoC = XC(end,:)'; 
    sols(c1).woutC=xwoC;       % equilibrium without viruses
    sols(c1).wout=xwoC./convC;
    % how far is the system from steady state at xwoC
    nrms(c1)=norm(virusremoval_model_rhsC(xwoC,pn));
end

%% Construct flux matrix and solve fluxes

% which parameter sets reached equilibrium
idxs = find(nrms < 1e-6);
nsel = numel(idxs);
psel = pars(idxs);

% preallocate only selected cases
flxref = struct('with', [], 'wout', []);
flxs = repmat(flxref, nsel, 1);
xeqref = struct('with', [], 'wout', []);
xeqsC = repmat(xeqref, nsel, 1);

for c2 = 1:nsel
    c1 = idxs(c2);
    pn = pars(c1);

    xeqC = sols(c1).withC;       % equilibrium with viruses
    [~, ph_with] = virusremoval_model_rhsC(xeqC, pn);
    ph_with = virusremoval_model_std(ph_with(:,1:7));

    xwoC = sols(c1).woutC;       % equilibrium without viruses
    [~, ph_wout] = virusremoval_model_rhsC(xwoC, pn);
    ph_wout = virusremoval_model_std(ph_wout(:,1:7));

    flxs(c2).with = ph_with;
    flxs(c2).wout = ph_wout;
    xeqsC(c2).with = xeqC;
    xeqsC(c2).wout = xwoC;
end

%% Save data
save('virusremoval_data.mat', 'flxs', 'psel', 'xeqsC');