function [SOL, SYS_PERS]=solve_adaptive_step(t, ...
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
% optimal time step has been calculated by
% the local controller during the last micro time step
dt_(end) = SYS_PERS.controller.h;

t_ = t(1);
tEnd = t(2);

rejected_iter_ = 0;
refined_iter_ = 0;
solver = SYS_PERS.sys.method_hdl;
% initialization of the persistant variables C_PERS_
% that contains the history of step size controller parameters
% of the taken micro timesteps (only the last macro time step is saved).
if( sysIndex_ == 1)
    
    C_PERS_.eEstVec = [];
    C_PERS_.rhofac = 0;
    
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
    
    % find how many micro timesteps are rejected
    i = length(C_PERS_.t);
    
    % find closest to the tEnd
    while tEnd < C_PERS_.t(i)
        i = i -1;
    end
    
    merge = 0;
    if (abs(tEnd - C_PERS_.t(i)) < 1e-7)
        i = i -1;
        merge  = 1;
    end
    
    % number of rejected micro timesteps
    rejected_iter_ = length(C_PERS_.t)-i;
    %number of total rejected micro timesteps
    
    sysIndex_ =  sysIndex_ - rejected_iter_;
    
    t_ = C_PERS_.t(i);
    %remove rejected solutions
    dt_ = dt_(1:end-rejected_iter_);
    y_ = y_(:,1:end-rejected_iter_);
    
    if (merge)
        dt_(end) = (tEnd - t_);
    end
    
    
    % restore the values from the history of the microtime steps
    C_PERS_.eEstVec = C_PERS_.eEstVec_(i,:);
    C_PERS_.rhofac = C_PERS_.rhofac_(i);
    
    C_PERS_.t = C_PERS_.t(1:i);
    C_PERS_.eEstVec_ = C_PERS_.eEstVec_(1:i,:);
    C_PERS_.rhofac_ = C_PERS_.rhofac_(1:i);
    
    
    S_PERS_.Fac = S_PERS_.fac_(:,i);
    S_PERS_.Delta_old = S_PERS_.delta_old_(:,:,i);
    S_PERS_.New_Jac = S_PERS_.new_jac_(i);
    S_PERS_.J = S_PERS_.j_(:,:,i);

    S_PERS_.fac_ = S_PERS_.fac_(:,1:i);
    S_PERS_.delta_old_ = S_PERS_.delta_old_(:,:,1:i);
    S_PERS_.new_jac_ = S_PERS_.new_jac_(1:i);
    S_PERS_.j_= S_PERS_.j_(:,:,1:i);
    
else
    % keep the history of micro timesteps only
    % from the last macro timestep
    C_PERS_.t = t_;
    C_PERS_.eEstVec_ = [];
    C_PERS_.eEstVec_ = C_PERS_.eEstVec;
    C_PERS_.rhofac_ = C_PERS_.rhofac;
    
    S_PERS_.fac_ = [];
    S_PERS_.fac_ = S_PERS_.Fac;
    S_PERS_.delta_old_ = [];
    S_PERS_.delta_old_ = S_PERS_.Delta_old;
    S_PERS_.new_jac_ = S_PERS_.New_Jac;
    S_PERS_.j_ = [];
    S_PERS_.j_= S_PERS_.J;
end

stat_=zeros(1,3);
stat=zeros(1,3);


% micro timestepping
rejected_ = false;

C_PERS_.eEst = -1;

while t_ < tEnd
    
    % Ensure that we dont leave the interval
    if (t_ + dt_(end) > tEnd)
        dt_(end) = tEnd-t_;
    end

    varsTilde = approximate(t_vars, y_vars,  -(t_-t(1)+dt_(end)));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), S_PERS_] = solver (...
        [t_ t_+dt_(end)],...
        y_, dt_,...
        SYS_PERS.sys.ode_hdl,...
        varsTilde, relTol, S_PERS_, rejected_);
    stat=stat+stat_;

    % calculate the error
    [eEst, ~] = ee_skelboe2000(...
        [y_ SOL], dt_, relTol,  SYS_PERS.solver.yTypical);
    
    % Update time and timestep
    [new_dt_,  rejected_, C_PERS_] = ec_h211b(...
        dt_, eEst, C_PERS_ );
    
    if (rejected_)
        dt_(end) = new_dt_;
        refined_iter_ = refined_iter_ + 1;
    else
        C_PERS_.eEst = max(C_PERS_.eEst, eEst);
        % TODO: allocatoin of memory should be considered
        t_ = t_ + dt_(end);
        % update controller history
        C_PERS_.t(end+1) = t_;
        C_PERS_.eEstVec_(end+1,:)= C_PERS_.eEstVec;
        C_PERS_.rhofac_(end+1) = C_PERS_.rhofac;
        % update solver history

        S_PERS_.fac_(:,end+1) = S_PERS_.Fac;
        S_PERS_.delta_old_(:,:,end+1) = S_PERS_.Delta_old;
        S_PERS_.new_jac_(end+1) = S_PERS_.New_Jac;
        S_PERS_.j_(:,:,end+1) = S_PERS_.J;
        % save solution
        dt_(end+1) = new_dt_;
        y_(:,end+1) = SOL;
        sysIndex_ = sysIndex_ + 1;
    end
    
end

% save the last calculated values for GS organization
SYS_PERS.sys.varsTilde = varsTilde;
% save solution
SYS_PERS.sol.y(:,1:sysIndex_) = y_;
SYS_PERS.sol.dt(1:sysIndex_-1) = dt_(1:end-1);
% save controller workspace
SYS_PERS.controller = C_PERS_;
SYS_PERS.controller.h = new_dt_;
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