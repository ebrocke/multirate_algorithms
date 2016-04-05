% <v> <t> vector values are filled in the following order: 
% from current time t(1) = 0, t(2) = h_{n}, t(3)=h_{n}+h_{n-1} and etc,
% where h_{n} is the last step size taken by a system.
function [v t SYSTEM] = get_exchanged_vector (GET_EXCH_HDL, STEP_REJECTED, SYSTEM, PARAMS)
global MODE

sysIndex = SYSTEM.stats.acceptedIter;

% if the last step has been rejected
% use the last exchanged vector
if STEP_REJECTED
    s_ = min(sysIndex,MODE);
    v = SYSTEM.sys.y_exch(1:s_,:);
    t = SYSTEM.sys.t_exch(1:s_);
    return;
end
% the size of the history (maximum 3)
s_ = min(sysIndex+1,MODE);
% retrieve system solution vector
y_ = SYSTEM.sol.y(:,1:sysIndex+1);
dt_= SYSTEM.sol.dt(1:sysIndex);
% number of exchanged variables coming from the system
N=2;
% allocate exchanged vector
v = zeros(s_,N);
t = [0 cumsum(fliplr(dt_(end-s_+2:end)))];
%fill in exchanged vector
jj_ = s_-1;
while jj_ >=0
    v_ = feval(GET_EXCH_HDL, SYSTEM, y_(:,end-jj_), t(jj_+1), PARAMS);
    v(jj_+1,:)    = v_;
    jj_ = jj_ -1;
end
% save exchanged vector
SYSTEM.sys.y_exch(1:s_,:) = v;
SYSTEM.sys.t_exch(1:s_) = t;

end