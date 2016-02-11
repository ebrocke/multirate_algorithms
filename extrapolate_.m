%% Extrapolation 0th,1st,2d order Skelboe(2000)
% y_tilde is matrix of maximum dimension three
% the last row 
function y_tilde = extrapolate(y, h)
global MODE

i_ = size(y,1);
y_tilde = zeros(min([i_+1,MODE]),size(y,2));
y_tilde(1:i_,:) = y;

if i_ == 1 % if first few values, use constant extrp.
    
    y_tilde(end,:) = y; 
    
elseif i_ == 2 % if later, use second order polyn.
    gamma = h(end)/h(end-1);  % = h(n)/h(n-1)
    
    y_tilde(end,:) = y(2,:)+gamma*(y(2,:)-y(1,:));
    
elseif i_ == 3
    y_tilde = circshift(y_tilde,-1);
    gamma = h(end)/h(end-1);
    delta = 1+h(end-2)/h(end-1);    % 1 + h(n-2)/h(n-1)
    a2 = gamma*(gamma+delta)/(1-delta);
    a3 = gamma*(gamma+1)/(delta*(delta-1));
    a1 = 1-a2-a3;
    
    y_tilde(end,:) = a1*y(3,:)+a2*y(2,:)+a3*y(1,:);
    
end

end