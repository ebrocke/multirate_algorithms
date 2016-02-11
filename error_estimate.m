function [ERR, ERR_INDEX] = error_estimate(Y, DT, SYS_INDEX, RELTOL, YTYPICAL)

in = Y (:,end); % the last element is the solution to evaluate

if SYS_INDEX > 2 % for later steps, use BDF2 with mode 3
    
    gamma = DT(end)/DT(end-1);
    delta = 1+DT(end-2)/DT(end-1);    % 1 + h(n-2)/h(n-1)
    a2 = gamma*(gamma+delta)/(1-delta);
    a3 = gamma*(gamma+1)/(delta*(delta-1));
    a1 = 1-a2-a3;
    Yp2n = a1*Y(:,end-1)...
        +a2*Y(:,end-2)...
        +a3*Y(:,end-3);
    
    [ERR, ERR_INDEX] = max(abs(in-Yp2n)./(RELTOL*(abs(in)+abs(YTYPICAL))));
    
elseif SYS_INDEX == 2 % For the second step, use BDF2 with mode 2
    
    gamma = DT(end)/DT(end-1);  % = h(n)/h(n-1)
    Ypn = Y(:,end-1)...
        +gamma*(Y(:,end-1)...
        -Y(:,end-2));
    
    [ERR, ERR_INDEX] = max(abs(in-Ypn)./(RELTOL*(abs(in)+abs(YTYPICAL))));
    
else
    ERR = 0;
    ERR_INDEX = -1;
end
end