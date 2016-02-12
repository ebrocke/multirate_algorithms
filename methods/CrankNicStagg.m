function [SOL, NJ_CALLS, NJ_ODE_CALLS, ODE_CALLS, PERSISTENT]= CrankNicStagg(T,...
    Y, ~, ~, ODE_PARAMS, ~, PERSISTENT, ~)
% Consider a function on the form
% x' = A(y)x + b(y)
% y' = c(x) + D(x)y

SOL = Y(:,end);
%[frac, cai] = ODE_PARAMS{:};

NJ_CALLS = 1;
NJ_ODE_CALLS = 0;
ODE_CALLS = 1;

% Defining constants.
T0 = T(1); % Starting time
TF = T(2); % Final time
h = (TF-T0); %/Nhh;

t = T0;
% only the knowledge of variables at t{n} are required
% to update the state variables at t{n+1}
% (constant interpolation over the interval is used)

t = t + h;

x = get_xvec(SOL); %gate probabilities
y = get_yvec(SOL); %voltages

% Update x
[A b] = getAb(SOL);
Id = eye(size(A));
if t == h % If first step, extrapolate the value at half step.
    x = x + A*x*h/2+b*h/2;
else % If not first, just updateS
    x = (Id-h*A/2)\((Id+h*A/2)*x + h*b);
    %         x = (Id-h*A)\((Id+h*A)*x + h*b);
end

SOL = set_xvec(SOL,x);

%Update y
[c D] = getcD(T, SOL, ODE_PARAMS{:});
Id = eye(size(D));
y = (Id-h*D/2)\((Id+h*D/2)*y + h*c);

SOL = set_yvec(SOL, y);

end

% function var_tilda = interpolate(vals)
% global dt
% if size(vals,1) == 1
%     var_tilda = vals(1);
% elseif size(vals,1)==2
%     %gamma = (t-t0)/dt;
%     var_tilda = (vals(2,1)+vals(1,1))/2;
% elseif size(vals,1)==3
%     gamma = 1;
%     delta = 3;
%     a2 = gamma*(gamma+delta)/(1-delta);
%     a3 = gamma*(gamma+1)/(delta*(delta-1));
%     a1 = 1-a2-a3;
%     var_tilda = (vals(3,1)-a2*vals(2,1)-a3*vals(1,1))/a1;
% else
%     var_tilda = 0;
%     disp('Something is wrong!')
% end
% end
