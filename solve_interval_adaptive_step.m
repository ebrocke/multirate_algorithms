function [SOL, SYS_PERS]=solve_interval_adaptive_step(solver, ...
    t, ...
    ode_fun, ...
    vars, ...
    relTol, ...
    SYS_PERS, ...
    step_rejected)

persistent SOLVER_PERS_
varsTilde = vars{1};
dtTilde = vars{2};

if (any(isnan(varsTilde)))
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    return;
end

sysIndex = SYS_PERS.stats.acceptedIter+1;
dt_ = SYS_PERS.sol.dt(1:sysIndex-1);
dt_(end+1) = SYS_PERS.controller.optimal_dt;
y_ = SYS_PERS.sol.y(:,1:sysIndex);

t_ = t(1);
tEnd = t(2);

% initialization of the persistant variables SOLVER_PERS_
% that contains the history of step size controller parameters
% of the taken micro timesteps (only the last macro time step is saved).
if( sysIndex == 1)
    
    SOLVER_PERS_.eEstVec = [NaN NaN NaN];
    SOLVER_PERS_.rhofac = 0;
    
    SOLVER_PERS_.eEstVec_= [];
    SOLVER_PERS_.rhofac_ = [];
end

% if the last macro timestep has been rejected
% delete the history up to the end of
% a new macro timestep, restore the state
% of the variables corresponding to this point in time.
if(step_rejected)
    
    % find how many micro timesteps are rejected
    i = length(SOLVER_PERS_.t);
    
    % find closest to the tEnd
    while tEnd < SOLVER_PERS_.t(i)
        i = i -1;
    end
    
    merge = 0;
    if (abs(tEnd - SOLVER_PERS_.t(i)) < 1e-10)
        i = i -1;
        merge  = 1;
    end
    
    % number of rejected micro timesteps
    k = length(SOLVER_PERS_.t)-i;
    %number of total rejected micro timesteps
    SYS_PERS.stats.rejectedIter =  SYS_PERS.stats.rejectedIter + k;
    SYS_PERS.stats.acceptedIter =  SYS_PERS.stats.acceptedIter - k;
    
    t_ = SOLVER_PERS_.t(i);
    %remove rejected solutions
    dt_ = dt_(1:end-k);
    
    if (merge)
        dt_(end) = (tEnd - t_);
    end
    
    y_ = y_(:,1:end-k);
    
    %t_ = SOLVER_PERS_.t(i);
    
    % update the history of the microtime steps
    SOLVER_PERS_.eEstVec = SOLVER_PERS_.eEstVec_(i,:);
    SOLVER_PERS_.rhofac = SOLVER_PERS_.rhofac_(i);
    
    SOLVER_PERS_.t = SOLVER_PERS_.t(1:i);
    SOLVER_PERS_.eEstVec_ = SOLVER_PERS_.eEstVec_(1:i,:);
    SOLVER_PERS_.rhofac_ = SOLVER_PERS_.rhofac_(1:i);
    
else
    % keep the history of micro timesteps only
    % from the last macro timestep
    SOLVER_PERS_.t = t_;
    SOLVER_PERS_.eEstVec_=[];
    SOLVER_PERS_.eEstVec_(1,:)= SOLVER_PERS_.eEstVec;
    SOLVER_PERS_.rhofac_ = SOLVER_PERS_.rhofac;
end

stat_=zeros(1,3);
stat=zeros(1,3);

% micro timestepping
rejected_ = false;

while t_ < tEnd
    
    % Ensure that we dont leave the interval
    if (t_ + dt_(end) > tEnd)
        dt_(end) = tEnd-t_;
    end

    varsTilde_ = interpolate(varsTilde, dtTilde, t_ + dt_(end) - t(1));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), SYS_PERS.solver] = solver ([t_ t_+dt_(end)],...
        y_, dt_,...
        ode_fun, varsTilde_, ...
        relTol, SYS_PERS.solver, rejected_);
    stat=stat+stat_;
    
    % Update time and timestep
    [new_dt,  rejected_, SOLVER_PERS_] = error_control([y_ SOL],...
        dt_, sysIndex, relTol, SYS_PERS.solver.yTypical, SOLVER_PERS_, 0 );
    
    if (rejected_)
        dt_(end) = new_dt;
        sysIndex = sysIndex - 1;
        SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + 1;
    else
        SYS_PERS.stats.acceptedIter = SYS_PERS.stats.acceptedIter + 1;
        t_ = t_ + dt_(end);
        SOLVER_PERS_.t(end+1) = t_;
        SOLVER_PERS_.eEstVec_(end+1,:)= SOLVER_PERS_.eEstVec;
        SOLVER_PERS_.rhofac_(end+1) = SOLVER_PERS_.rhofac;
        dt_(end+1) = new_dt;
        y_(:,end+1) = SOL;
    end
    sysIndex = sysIndex + 1;
end

sysIndex = SYS_PERS.stats.acceptedIter+1;
if(sysIndex ~= size(y_,2))
    stop = 0;
end
SYS_PERS.sol.y(:,1:sysIndex) = y_;
SYS_PERS.sol.dt(1:sysIndex-1) = dt_(1:end-1);
SYS_PERS.controller.optimal_dt = new_dt;
SYS_PERS.controller.err = SOLVER_PERS_.err;

SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);


end