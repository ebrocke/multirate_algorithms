function [SOL, SYS_PERS]=solve_interval_one_step(solver, ...
    t, ...
    ode_fun, ...
    vars, ...
    relTol, ...
    SYS_PERS, ...
    step_rejected)

H = t(2)-t(1);
if (step_rejected)
   SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + 1;   
else
   SYS_PERS.stats.acceptedIter = SYS_PERS.stats.acceptedIter + 1;
end
varsTilde = vars{1};
if (any(isnan(varsTilde)))
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    return;
end
sysIndex = SYS_PERS.stats.acceptedIter;

 

SYS_PERS.sol.dt(sysIndex) = H; % optimal predicted time step

stat=zeros(1,3);

[SOL, stat(1), stat(2), stat(3), SYS_PERS.solver] = solver(t,...
    SYS_PERS.sol.y(:,1:sysIndex), ...
    SYS_PERS.sol.dt(1:sysIndex), ...
    ode_fun, varsTilde(end,:), relTol, SYS_PERS.solver, step_rejected);

SYS_PERS.sol.y(:, sysIndex+1) = SOL; % save the solution

SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);
end