function [SOL, SYS_PERS]=solve_one_step(t, ...
    vars, ...
    relTol, ...
    step_rejected, ...
    SYS_PERS)

% update statistics
if (step_rejected)
   SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + 1;
else
   SYS_PERS.stats.acceptedIter = SYS_PERS.stats.acceptedIter + 1;   
end

% return NaN if previous solver failed
varsTilde = vars{2};
if (any(isnan(varsTilde)))
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    return;
end

stat=zeros(1,3);

sysIndex = SYS_PERS.stats.acceptedIter;

solver = SYS_PERS.sys.method_hdl;

SYS_PERS.sol.dt(sysIndex) = t(2)-t(1); % optimal predicted time step

[SOL, stat(1), stat(2), stat(3), SYS_PERS.solver] = solver(t,...
    SYS_PERS.sol.y(:,1:sysIndex), ...
    SYS_PERS.sol.dt(1:sysIndex), ...
    SYS_PERS.sys.ode_hdl,...
    varsTilde(1,:), relTol, SYS_PERS.solver, step_rejected);

SYS_PERS.sol.y(:, sysIndex+1) = SOL; % save the solution


% update statistics
SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);
end