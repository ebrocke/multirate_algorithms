function [H_MICRO, PERSISTENT] = ec_classical(H_MACRO, E_EST, PERSISTENT)

hmax = 1e-2;%1e-5;
mmax = 2;
p = 1./2.;

if (PERSISTENT.init > 0)
    
    m_ = PERSISTENT.m;
    h_ = PERSISTENT.h;
    
    if (E_EST < 1)
        h_ = min(hmax, h_ / max(0.2, 1.25*E_EST^(p)));
        m_ = max(fix(m_/1.5),...
            max(mmax,fix(H_MACRO/h_)));
    else
        h_ = h_*max(.1,.5*E_EST^(-p));
        m_ = min(fix(1.5*m_),...
            max(mmax,fix(H_MACRO/h_)));
    end
else
    PERSISTENT.init = 1;
    m_ = mmax;
end
H_MICRO=H_MACRO/m_;
PERSISTENT.m = m_;
PERSISTENT.h = H_MICRO;
end