function [SOL, SYS_PERS]=solve_fixed_step(t, ...
    vars, ...
    relTol, ...
    step_rejected,...
    SYS_PERS)

t_vars = vars{1};
y_vars = vars{2};

if (any(isnan(y_vars)))
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    return;
end

% controller workspace
C_PERS_ = SYS_PERS.controller;
% solver workspace
S_PERS_ = SYS_PERS.solver;

sysIndex_ = SYS_PERS.stats.acceptedIter+1;

% solution vector
y_ = SYS_PERS.sol.y(:,1:sysIndex_);
dt_ = SYS_PERS.sol.dt(1:sysIndex_);


t_ = t(1);
tEnd = t(2);


rejected_iter_ = 0;
refined_iter_ = 0;
solver = SYS_PERS.sys.method_hdl;

if (sysIndex_ == 1)
    S_PERS_.Fac = [];
    S_PERS_.Delta_old = [];
    S_PERS_.New_Jac = true;
    S_PERS_.J = [];
end

% if the last macro timestep has been rejected
% delete the history up to the end of
% a new macro timestep, restore the state
% of the variables corresponding to this point in time.
if(step_rejected)
 
    %number of total rejected micro timesteps  
    rejected_iter_ = C_PERS_.m;
    
    sysIndex_ =  sysIndex_ - rejected_iter_;
   
    %remove rejected solutions
    dt_ = dt_(1:end-rejected_iter_);
    y_ = y_(:,1:end-rejected_iter_);
    
    % resotore solver workspace
    S_PERS_.Fac = S_PERS_.fac_;
    S_PERS_.Delta_old = S_PERS_.delta_old_;
    S_PERS_.New_Jac = S_PERS_.new_jac_;
    S_PERS_.J = S_PERS_.j_;
    
else
    % save the solver workspace
    S_PERS_.fac_ = S_PERS_.Fac;
    S_PERS_.delta_old_ = S_PERS_.Delta_old;
    S_PERS_.new_jac_ = S_PERS_.New_Jac;
    S_PERS_.j_= S_PERS_.J;
 
end

stat_=zeros(1,3);
stat=zeros(1,3);

% suggested by the global controller
dt_(end) = SYS_PERS.sys.h;

% micro timestepping
rejected_ = false;
jj_ = SYS_PERS.sys.m;

while jj_ > 0
    
    varsTilde = approximate(t_vars, y_vars,  -(t_-t(1)+dt_(end)));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), S_PERS_] = solver (...
        [t_ t_+dt_(end)],...
        y_, dt_,...
        SYS_PERS.sys.ode_hdl,...
        varsTilde, relTol, S_PERS_, rejected_);
    stat=stat+stat_;
    %dbstop if warning
    % allocatoin of memory should be considered
    t_ = t_ + SYS_PERS.sys.h;
    y_(:,end+1) = SOL;
    dt_(end+1) = SYS_PERS.sys.h;
    sysIndex_ = sysIndex_ + 1;
    jj_ = jj_ - 1;
    if any(isnan(SOL))
        break;
    end
    
end
C_PERS_.m = SYS_PERS.sys.m-jj_;
% save the last calculated values for GS organization
SYS_PERS.sys.varsTilde = varsTilde;
% calculate the error
[C_PERS_.eEst, ~] = ee_skelboe2000(...
        y_, dt_(1:end-1), relTol,  SYS_PERS.solver.yTypical);
% save solution
SYS_PERS.sol.y(:,1:sysIndex_) = y_;
SYS_PERS.sol.dt(1:sysIndex_-1) = dt_(1:end-1);
% save controller workspace
SYS_PERS.controller = C_PERS_;
% save solver workspace
SYS_PERS.solver = S_PERS_;
% update statistics
SYS_PERS.stats.acceptedIter = sysIndex_-1;
SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + refined_iter_;
SYS_PERS.stats.rejectedIter = SYS_PERS.stats.rejectedIter + rejected_iter_;
SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);

end
