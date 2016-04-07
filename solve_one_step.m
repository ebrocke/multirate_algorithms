function [SOL, SYS_PERS]=solve_one_step(t, ...
    vars, ...
    relTol, ...
    SYS_PERS)

% return NaN if previous solver failed
%varsTilde = vars{2};
solve = true;
if (any(isnan(vars{2})))
    %display('solver_one_step:: varsTilde is NaN')
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    solve = false;
   % return;
end

% aligns solution according to interval t
SYS_PERS = update(t, SYS_PERS);
% controller workspace
C_PERS_ = SYS_PERS.controller;
% solver workspace
S_PERS_ = SYS_PERS.solver;

stat=zeros(1,3);

sysIndex_ = SYS_PERS.stats.acceptedIter+1;

% solution vector
y_ = SYS_PERS.sol.y(:,1:sysIndex_);
dt_ = SYS_PERS.sol.dt(1:sysIndex_);

% optimal time step has been calculated by
% the local controller during the last micro time step
%dt_(end) = SYS_PERS.controller.h;

solver = SYS_PERS.sys.method_hdl;
t_ = C_PERS_.t(end);
h_ = dt_(end);

if solve
    varsTilde = approximate(vars{1}, vars{2}, -(t(2)-t(1)));
    [SOL, stat(1), stat(2), stat(3), S_PERS_] = solver([t_ t_+h_],...
        y_, dt_, ...
        SYS_PERS.sys.ode_hdl,...
        varsTilde(1,:), relTol, S_PERS_);
end

% calculate the error
[eEst, ~] = ee_skelboe2000(...
    [y_ SOL], dt_, relTol,  S_PERS_.yTypical);

% Update time and timestep
[h_,  S_PERS_.step_rejected, C_PERS_] = C_PERS_.fn(...
    dt_, eEst, C_PERS_);

% update statistics
if (S_PERS_.step_rejected)
   S_PERS_.step_rejected = false; % we restore the values outside
   SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + 1;
   SOL = NaN(size(SYS_PERS.solver.yTypical));
   SYS_PERS.sol.dt(sysIndex_) = h_;
else
   SYS_PERS.stats.acceptedIter = SYS_PERS.stats.acceptedIter + 1;
   SYS_PERS.sol.y(:, sysIndex_+1) = SOL; % save the solution 
   SYS_PERS.sol.dt(sysIndex_+1) = h_;
   %SYS_PERS.sol.dt(sysIndex_) = dt_(end);
   
   % update controller history
   C_PERS_.t(end+1) = t_+dt_(end);
   C_PERS_.eEstVec_(end+1,:)= C_PERS_.eEstVec;
   C_PERS_.rhofac_(end+1) = C_PERS_.rhofac;
   C_PERS_.eEst_(end+1) = eEst;
   % update solver history
   S_PERS_.fac_(:,end+1) = S_PERS_.Fac;
   S_PERS_.delta_old_(:,:,end+1) = S_PERS_.Delta_old;
   S_PERS_.new_jac_(end+1) = S_PERS_.New_Jac;
   S_PERS_.j_(:,:,end+1) = S_PERS_.J;
end

SYS_PERS.controller = C_PERS_;
SYS_PERS.solver = S_PERS_;
%SYS_PERS.controller.h = h_;
SYS_PERS.controller.eEst = eEst;
% update statistics
SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);
end