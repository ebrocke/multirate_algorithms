% save the history of controller and solver parameters
% only for the current interval t (macrotime step)
function [SYSTEM] = update(t, SYSTEM)
tEnd = t(2);
% controller workspace
C_PERS_ = SYSTEM.controller;
% solver workspace
S_PERS_ = SYSTEM.solver;

sysIndex_ = SYSTEM.stats.acceptedIter+1;
if sysIndex_ == 1
    return
end

% solution vector
y_ = SYSTEM.sol.y(:,1:sysIndex_);
dt_ = SYSTEM.sol.dt(1:sysIndex_);
%dt_(end) = SYSTEM.controller.h;

ii_ = length(C_PERS_.t);
% if (ii_ > 2)
%     stop = 0;
% end
% find index closest to the end time
%(on the left from the end)
while tEnd < C_PERS_.t(ii_)
    ii_ = ii_ -1;
end

% % we merge intervals if the step size less than 1e-7
% % a good check if something wrong in the integration
% % since the step size should not go below 1e-7
% merge = 0;
% if (abs(tEnd - C_PERS_.t(ii_)) < 1e-7)
%     ii_ = ii_ -1;
%     merge  = 1;
% end

% number of rejected micro timesteps
rejected_iter_ = length(C_PERS_.t)-ii_;

sysIndex_ =  sysIndex_ - rejected_iter_;

% current time
%t_ = C_PERS_.t(ii_);

%remove rejected solutions
dt_ = dt_(1:end-rejected_iter_);
y_ = y_(:,1:end-rejected_iter_);

% adjust optimal time step
% if (merge)
%     dt_(end) = (tEnd - t_);
% end

% find index closest to the start time
%(on the left from the start)
jj_ =  length(C_PERS_.t)-rejected_iter_;
while jj_>1 && t(1) < C_PERS_.t(jj_)
    jj_ = jj_ -1;
end

% we save the history only for the current macro time step
% if the last macro step rejected shrink the history
% otherwise reassign with the last values
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

% save solution
SYSTEM.sol.y(:,1:sysIndex_) = y_;
SYSTEM.sol.dt(1:sysIndex_) = dt_(1:end);
% save controller workspace
SYSTEM.controller = C_PERS_;
%SYSTEM.controller.h = dt_(end);
% save solver workspace
SYSTEM.solver = S_PERS_;
% update statistics
SYSTEM.stats.acceptedIter = sysIndex_-1;
SYSTEM.stats.rejectedIter = SYSTEM.stats.rejectedIter + rejected_iter_;
end