function [H_MICRO, PERSISTENT] = ec_h211b_hmicro(H_MACRO, E_EST, PERSISTENT)

mmax = 2;

if (PERSISTENT.init == 0)
    m_ = mmax;
    PERSISTENT.h = zeros(1,3);
    PERSISTENT.h(end) = 1e-5;
else
    m_ = PERSISTENT.m;
end
dt_ = PERSISTENT.h;
[h_,  step_rejected_, PERSISTENT] = ec_h211b(dt_, E_EST, PERSISTENT);
if (step_rejected_)
    m_ = min(fix(1.15*m_),...
        max(mmax,fix(H_MACRO/h_)));
  
else   
    m_ = max(fix(m_/1.5),...
        max(mmax,fix(H_MACRO/h_)));
end

H_MICRO=H_MACRO/m_;
PERSISTENT.m = m_;
PERSISTENT.h = circshift(dt_,[0,-1]);
PERSISTENT.h(end) = H_MICRO;
end