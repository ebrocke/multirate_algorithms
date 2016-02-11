%% Point interpolation ( Skelboe(2000) representation)
% y(,:) is a vector of solutions <y(1,:) y(2,:) y(3,:)>
% obtained with time steps <h(1) h(2)> respectively
%<y_tilde> interpolated solution at dt \in h(2)

%         h(1)      dt     h(2)
%   *------------*------*---------------*
%   y(1,:)      y(2,:)  y_tilde       y(3,:)

function y_tilde = interpolate(dt, y, h)

y_tilde = zeros(1,size(y,2));

if size(y,1) == 1 % if first few values, use constant extrp.
    
    y_tilde = y;
    
elseif size(y,1) == 2 % if later, 1st order interpolation
    gamma = (h(end)-dt)/dt;  % = h(n)/h(n-1)
    
    y_tilde(:) = (y(2,:)+gamma*y(1,:))/(1+gamma);
    
elseif size(y,1) == 3
    gamma = (h(end)-dt)/dt;
    delta = 1+h(end-1)/dt;    % 1 + h(n-2)/h(n-1)
    a2 = gamma*(gamma+delta)/(1-delta);
    a3 = gamma*(gamma+1)/(delta*(delta-1));
    a1 = 1-a2-a3;
    
    y_tilde(:) = (y(3,:)-a2*y(2,:)-a3*y(1,:))/a1;
    
end
end
