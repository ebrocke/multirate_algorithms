function [SOL, SYS_PERS]=solve_adaptive_step(t, ...
    vars, ...
    relTol, ...
    SYS_PERS)

t_vars = vars{1};
y_vars = vars{2};

if (any(isnan(y_vars)))
    SOL = NaN(size(SYS_PERS.sol.y,1),1);
    return;
end
% aligns solutions according to interval t
% updates the history
SYS_PERS = update(t, SYS_PERS);

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
%dt_(end) = SYS_PERS.controller.h;

tEnd = t(2);
refined_iter_ = 0;
solver = SYS_PERS.sys.method_hdl;

stat_=zeros(1,3);
stat=zeros(1,3);

% last value is an optimal time step
h_ = dt_(end);
% current time
t_ = C_PERS_.t(end);

while t_ < tEnd
 
    %% use for synch comm, comment for usynch
    % Ensure that we dont leave the interval
    if (t_ + h_ > tEnd)
        dt_(end) = tEnd-t_;
    end
    %% end
    
    varsTilde = approximate(t_vars, y_vars,  -(t_-t(1)+h_));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), S_PERS_] = solver (...
        [t_ t_+h_],...
        y_, dt_,...
        SYS_PERS.sys.ode_hdl,...
        varsTilde, relTol, S_PERS_);
    stat=stat+stat_;
    
    % calculate the error
    [eEst, ~] = ee_skelboe2000(...
        [y_ SOL], dt_, relTol,  SYS_PERS.solver.yTypical);
    
    % Update time and timestep
    [h_,  S_PERS_.step_rejected, C_PERS_] = C_PERS_.fn(...
        dt_, eEst, C_PERS_ );
    
    if (S_PERS_.step_rejected)
        dt_(end) = h_;
        refined_iter_ = refined_iter_ + 1;
    else
        % TODO: allocatoin of memory should be considered
        t_ = t_ + dt_(end);
        % update controller history
        C_PERS_.t(end+1) = t_;
        C_PERS_.eEstVec_(end+1,:)= C_PERS_.eEstVec;
        C_PERS_.rhofac_(end+1) = C_PERS_.rhofac;
        C_PERS_.eEst_(end+1) = eEst;
        % update solver history
        S_PERS_.fac_(:,end+1) = S_PERS_.Fac;
        S_PERS_.delta_old_(:,:,end+1) = S_PERS_.Delta_old;
        S_PERS_.new_jac_(end+1) = S_PERS_.New_Jac;
        S_PERS_.j_(:,:,end+1) = S_PERS_.J;
        % save solution
        dt_(end+1) = h_;
        y_(:,end+1) = SOL;
        sysIndex_ = sysIndex_ + 1;
    end
    
end
%use for usynchronized comm, comment for synch
% if (t_ > tEnd)
%     delta_ = C_PERS_.t(end)-tEnd;
%     ts = [0 dt_(end-1) dt_(end-1)+dt_(end-2)];
%     ys = fliplr(y_(:,end-2:end));
%     SOL = approximate(ts,ys',delta_);
%     SOL = SOL';
%     %varsTilde = approximate(t_vars, y_vars,  -(t(2)-t(1)));
% end
%
% save the last calculated values for GS organization
%SYS_PERS.sys.varsTilde = varsTilde;
% save solution
SYS_PERS.sol.y(:,1:sysIndex_) = y_;
SYS_PERS.sol.dt(1:sysIndex_) = dt_(1:end);
% save controller workspace
SYS_PERS.controller = C_PERS_;
%SYS_PERS.controller.h = h_;
SYS_PERS.controller.eEst = max(C_PERS_.eEst_);
% save solver workspace
SYS_PERS.solver = S_PERS_;
% update statistics
SYS_PERS.stats.acceptedIter = sysIndex_-1;
SYS_PERS.stats.refinedIter = SYS_PERS.stats.refinedIter + refined_iter_;
SYS_PERS.stats.numjac = SYS_PERS.stats.numjac + stat(1);
SYS_PERS.stats.n_ode_numjac = SYS_PERS.stats.n_ode_numjac + stat(2);
SYS_PERS.stats.n_ode_iter = SYS_PERS.stats.n_ode_iter + stat(3);

end