% <v> <t> vector values are filled in the following order: 
% from current time t(1) = 0, t(2) = h_{n}, t(3)=h_{n}+h_{n-1} and etc,
% where h_{n} is the last step size taken by a system.
function [v t SYSTEM] = get_exchanged_vector (SYSTEM, STEP_REJECTED, GET_EXCH_HDL)
global MODE

sysIndex = SYSTEM.stats.acceptedIter+1;
y_ = SYSTEM.sol.y(:,1:sysIndex);
dt_= SYSTEM.sol.dt(1:sysIndex-1);

% the size of the history (maximum 3)
s_ = min(sysIndex,MODE);

% if the last step has been rejected
% use the last exchanged vector
if STEP_REJECTED
    v = SYSTEM.sys.y_exch(sysIndex - 1:end,:);
    t = SYSTEM.sys.dt_exch(sysIndex - 1:end);
    return;
end

% number of exchanged variables coming from the system
N=2;
% allocate exchanged vector
v = zeros(s_,N);
t = [0 cumsum(fliplr(dt_(end-s_+2:end)))];
%fill in exchanged vector
jj_ = s_-1;
while jj_ >=0
    v_ = feval(GET_EXCH_HDL, y_(:,end-jj_), t(jj_+1));
    v(jj_+1,:)    = v_;
    jj_ = jj_ -1;
end
% save exchanged vector
SYSTEM.sys.y_exch(end-s_+1:end,:) = v;
SYSTEM.sys.dt_exch(end-s_+1:end) = t;


end