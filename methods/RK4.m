function [SOL, NJ_CALLS, NJ_ODE_CALLS, ODE_CALLS, PERSISTENT] = RK4(T,...
    Y, ~, ODE_FUN, ODE_PARAMS, ~, PERSISTENT, ~)

NJ_CALLS = 0;
NJ_ODE_CALLS = 0;
ODE_CALLS = 4;

%Y = RK4(t, erk_vals, cell_vals, get_exch_vals, pproc_hdl,erk_ext)
%global Nerk
h = (T(2)-T(1)); %/(Nerk);
% Same outer step size as inner
[~, ca_flux] = ODE_PARAMS{:};
% ca_flux = pproc_hdl(t, ca_flux);
SOL = Y(:,end);

% options  = odeset('Reltol',1e-2,'Abstol',1e-9,...% 'MaxOrder', 2,...
%     'Stats','off');
t = T(1);
%for ii = 1:Nerk
k1 = ODE_FUN(t, SOL, 0, ca_flux);
t_ = t + h/2;
% ca_flux = extrapolate(t_, t(1), ca_flux);
k2 = ODE_FUN(t_, SOL+h*k1/2, 0, ca_flux);
k3 = ODE_FUN(t_, SOL+h*k2/2, 0, ca_flux);
t = t + h;
%ca_flux = extrapolate(t, t(1), ca_flux);
k4 = ODE_FUN(t, SOL+h*k3, 0, ca_flux);
SOL = SOL+h*(k1+2*k2+2*k3+k4)/6;
%end
%Y = Y(:,end);
end

% function var_tilda = extrapolate(t, t0, vals)
%     if size(vals,1) == 1
%         var_tilda = vals(1);
%     elseif size(vals,1)==2
%         gamma = (t-t0)/dt;
%         var_tilda = vals(2,1)+gamma*(vals(2,1)-vals(1,1));
%     elseif size(vals,1)==3
%         gamma = (t-t0)/dt;
%         delta = 2;
%         a2 = gamma*(gamma+delta)/(1-delta);
%         a3 = gamma*(gamma+1)/(delta*(delta-1));
%         a1 = 1-a2-a3;
%         var_tilda = a1*vals(3,1)+a2*vals(2,1)+a3*vals(1,1);
%     else
%         var_tilda = 0;
%         disp('Something is wrong!')
%     end
% end

