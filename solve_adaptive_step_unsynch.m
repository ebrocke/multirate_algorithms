function [SOL, SYS_PERS]=solve_adaptive_step_unsynch(t, ...
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

refined_iter_ = 0;
solver = SYS_PERS.sys.method_hdl;
% initialization of the persistant variables C_PERS_
% that contains the history of step size controller parameters
% of the taken micro timesteps (only the last macro time step is saved).
if( sysIndex_ == 1)
    
    C_PERS_.eEstVec = [NaN NaN NaN];
    C_PERS_.rhofac = 0;
    C_PERS_.t = t(1);
    C_PERS_.eEst_ = 0;
    
    S_PERS_.Fac = zeros(23,1);
    S_PERS_.Delta_old = zeros(23,3);
    S_PERS_.New_Jac = true;
    S_PERS_.J = zeros(23,23);
    
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


% find how many micro timesteps are rejected
ii_ = length(C_PERS_.t);

% find index closest to the end time
%(on the left from the end)
while t(2) < C_PERS_.t(ii_)
    ii_ = ii_ -1;
end

% number of rejected micro timesteps
rejected_iter_ = length(C_PERS_.t)-ii_;

sysIndex_ =  sysIndex_ - rejected_iter_;

% current time
t_ = C_PERS_.t(ii_);

% remove rejected solutions
dt_ = dt_(1:end-rejected_iter_);
y_ = y_(:,1:end-rejected_iter_);

% find index closest to the start time
%(on the left from the start)
jj_ =  length(C_PERS_.t)-rejected_iter_;
while t(1) < C_PERS_.t(jj_)
    jj_ = jj_ -1;
end

% we save the history only for the current integrtion interval
C_PERS_.eEstVec = C_PERS_.eEstVec_(ii_,:);
C_PERS_.rhofac = C_PERS_.rhofac_(ii_);

C_PERS_.t = C_PERS_.t(jj_:ii_);
C_PERS_.eEst_ = C_PERS_.eEst_(jj_:ii_);
C_PERS_.eEstVec_ = C_PERS_.eEstVec_(jj_:ii_,:);
C_PERS_.rhofac_ = C_PERS_.rhofac_(jj_:ii_);

S_PERS_.Fac = S_PERS_.fac_(:,ii_);
S_PERS_.Delta_old = S_PERS_.delta_old_(:,:,ii_);
S_PERS_.New_Jac = S_PERS_.new_jac_(ii_);
S_PERS_.J = S_PERS_.j_(:,:,ii_);

S_PERS_.fac_ = S_PERS_.fac_(:,jj_:ii_);
S_PERS_.delta_old_ = S_PERS_.delta_old_(:,:,jj_:ii_);
S_PERS_.new_jac_ = S_PERS_.new_jac_(jj_:ii_);
S_PERS_.j_= S_PERS_.j_(:,:,jj_:ii_);

% init
stat_=zeros(1,3);
stat=zeros(1,3);

rejected_ = false;
h_ = dt_(end);

while t_ < t(2)
    
    % extrapolation of slow variables
    varsTilde = approximate(t_vars, y_vars,  -(t_-t(1)+h_));
    
    % Calculate the solution
    [SOL, stat_(1), stat_(2), stat_(3), S_PERS_] = solver (...
        [t_ t_+h_],...
        y_, dt_,...
        SYS_PERS.sys.ode_hdl,...
        varsTilde, relTol, S_PERS_, rejected_);
    stat=stat+stat_;
    
    % calculate the error
    [eEst, ~] = ee_skelboe2000(...
        [y_ SOL], dt_, relTol,  SYS_PERS.solver.yTypical);
    
    %err = max(eEst, SYS_PERS.sys.eEst);
    % Update time and timestep
    [h_,  rejected_, C_PERS_] = ec_h211b(...
        dt_, eEst, C_PERS_ );
    
    if (rejected_)
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
        ii_ = ii_+1;
    end
    
end
if (t_ > t(2))
    delta_ = C_PERS_.t(end)-t(2);
    ts = [0 dt_(end-1) dt_(end-1)+dt_(end-2)];
    ys = fliplr(y_(:,end-2:end));
    SOL = approximate(ts,ys',delta_);
    SOL = SOL';
    %varsTilde = approximate(t_vars, y_vars,  -(t(2)-t(1)));
end

% save the last calculated values for GS organization
%SYS_PERS.sys.varsTilde = varsTilde;
% save solution
SYS_PERS.sol.y(:,1:sysIndex_) = y_;
SYS_PERS.sol.dt(1:sysIndex_-1) = dt_(1:end-1);
% save controller workspace
SYS_PERS.controller = C_PERS_;
SYS_PERS.controller.h = h_;

SYS_PERS.controller.eEst = max(C_PERS_.eEst_);
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