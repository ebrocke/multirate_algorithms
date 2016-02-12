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
SOLVER_PERS_ = SYS_PERS.controller;

sysIndex_ = SYS_PERS.stats.acceptedIter+1;

% solution vector
y_ = SYS_PERS.sol.y(:,1:sysIndex_);
dt_ = SYS_PERS.sol.dt(1:sysIndex_);


t_ = t(1);
tEnd = t(2);



rejected_iter_ = 0;
refined_iter_ = 0;
solver = SYS_PERS.sys.method_hdl;
% initialization of the persistant variables SOLVER_PERS_
% that contains the history of step size controller parameters
% of the taken micro timesteps (only the last macro time step is saved).
if( sysIndex_ == 1)
    
    SOLVER_PERS_.eEstVec = [NaN NaN NaN];
    SOLVER_PERS_.rhofac = 0;
    
    %SOLVER_PERS_.eEstVec_= [];
    %SOLVER_PERS_.rhofac_ = [];
end

% if the last macro timestep has been rejected
% delete the history up to the end of
% a new macro timestep, restore the state
% of the variables corresponding to this point in time.
if(step_rejected)
    
    %     % find how many micro timesteps are rejected
    %     i = length(SOLVER_PERS_.t);
    %
    %     % find closest to the tEnd
    %     while tEnd < SOLVER_PERS_.t(i)
    %         i = i -1;
    %     end
    %
    %     merge = 0;
    %     if (abs(tEnd - SOLVER_PERS_.t(i)) < 1e-10)
    %         i = i -1;
    %         merge  = 1;
    %     end
    %
    %     % number of rejected micro timesteps
    %     rejected_iter_ = length(SOLVER_PERS_.t)-i;
    %number of total rejected micro timesteps
    
    rejected_iter_ = SOLVER_PERS_.m;
    
    sysIndex_ =  sysIndex_ - rejected_iter_;
    
    %    t_ = SOLVER_PERS_.t(i);
    %remove rejected solutions
    dt_ = dt_(1:end-rejected_iter_);
    y_ = y_(:,1:end-rejected_iter_);
    
    %     if (merge)
    %         dt_(end) = (tEnd - t_);
    %     end
    
    
    % update the history of the microtime steps
    %SOLVER_PERS_.eEstVec = SOLVER_PERS_.eEstVec_(i,:);
    %SOLVER_PERS_.rhofac = SOLVER_PERS_.rhofac_(i);
    
    %SOLVER_PERS_.t = SOLVER_PERS_.t(1:i);
    %SOLVER_PERS_.eEstVec_ = SOLVER_PERS_.eEstVec_(1:i,:);
    %SOLVER_PERS_.rhofac_ = SOLVER_PERS_.rhofac_(1:i);
    
    %else
    % keep the history of micro timesteps only
    % from the last macro timestep
    %SOLVER_PERS_.t = t_;
    %SOLVER_PERS_.eEstVec_=[];
    %SOLVER_PERS_.eEstVec_(1,:)= SOLVER_PERS_.eEstVec;
    %SOLVER_PERS_.rhofac_ = SOLVER_PERS_.rhofac;
end

stat_=zeros(1,3);
stat=zeros(1,3);

% suggested by the global controller
jj_ = fix((tEnd-t_)/SYS_PERS.sys.h);
dt_(end) = SYS_PERS.sys.h;

% micro timestepping
rejected_ = false;
SOLVER_PERS_.m = jj_;
sysIndex_ = sysIndex_ + jj_;
while jj_ > 0
    
    
    varsTilde = approximate(t_vars, y_vars,  -(t_-t(1)+dt_(end)));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), SYS_PERS.solver] = solver (...
        [t_ t_+dt_(end)],...
        y_, dt_,...
        SYS_PERS.sys.ode_hdl,...
        varsTilde, relTol, SYS_PERS.solver, rejected_);
    stat=stat+stat_;
    
    % allocatoin of memory should be considered
    t_ = t_ + SYS_PERS.sys.h;
    y_(:,end+1) = SOL;
    dt_(end+1) = SYS_PERS.sys.h;
    
    jj_ = jj_ - 1;
end
% calculate the error
[SOLVER_PERS_.eEst, ~] = ee_skelboe2000(...
        [y_ SOL], dt_(1:end-1), sysIndex_-1, relTol,  SYS_PERS.solver.yTypical);
% save solution
SYS_PERS.sol.y(:,1:sysIndex_) = y_;
SYS_PERS.sol.dt(1:sysIndex_-1) = dt_(1:end-1);
% save controller workspace
SYS_PERS.controller = SOLVER_PERS_;
% update statistics
SYS_PERS.stats.acceptedIter = sysIndex_-1;
SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + refined_iter_;
SYS_PERS.stats.rejectedIter = SYS_PERS.stats.rejectedIter + rejected_iter_;
SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);

end
