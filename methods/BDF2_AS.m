function [SOL, NJ_CALLS, NJ_ODE_CALLS, ODE_CALLS, PERSISTENT] = BDF2_AS(T,...
    Y, DT, ODE_FUN, ODE_PARAMS, relTol, PERSISTENT, STEP_REJECTED)


Y_TYPICAL = PERSISTENT.yTypical;
if PERSISTENT.init
    FAC=[];
    DELTA_OLD = zeros(size(Y_TYPICAL,1), 3);
    PERSISTENT.init = false;
else
    FAC = PERSISTENT.Fac;
    DELTA_OLD = PERSISTENT.Delta_old;
end

if STEP_REJECTED
    %shift columns to the right to erase last solution
    DELTA_OLD = circshift(DELTA_OLD,[0,1]);
    %DELTA_OLD(:,1) = zeros(size(Y_TYPICAL));
end

NJ_CALLS = 0;
NJ_ODE_CALLS = 0;
ODE_CALLS = 0;

params = num2cell(ODE_PARAMS);

epsilon = 0.1*relTol*norm(Y_TYPICAL,Inf);
itmax = 5;
iter = 0;
h = T(2)-T(1);
t = T(2); % one step

if t == h || size(Y,2) == 2
    % Implicit Euler step
    Y_old = Y(:,end);
    
     g = @(t,Y_new) bsxfun(@minus,bsxfun(@minus,Y_new,Y_old),...
         DT(end)*feval(ODE_FUN,t,Y_new, params{:}));
    %g = @(t,Y_new) Y_new-Y_old-DT(end)*feval(ODE_FUN,t,Y_new, ODE_PARAMS{:});

    
    % Constant interpolation
    Y_guess = Y_old;
    
    % Prepare for the modified newton iteration
    gY = g(t,Y_guess);
    ODE_CALLS = ODE_CALLS + 1;
    
    [J,FAC,~,nf] = numjac(g,t,Y_guess,gY,{1e-6,Y_TYPICAL}, FAC,0);
    NJ_CALLS = NJ_CALLS + 1;
    NJ_ODE_CALLS = NJ_ODE_CALLS + nf;
    [LJ, UJ] = lu(J);
    
    % Solve via iteration
    while (iter < itmax)
        iter = iter + 1;
        
        dY = UJ\(LJ\gY);
                
        newnrm = norm(dY,Inf);
        if iter > 1 && newnrm >= 0.9*oldnrm
            %disp(['Slowly convergent Newton method, hh, h = ', num2str(dtCellVec(end))])
            SOL = NaN(size(Y_TYPICAL));
            DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
            %PERSISTENT=[{FAC} {DELTA_OLD}];
            PERSISTENT.Fac = FAC;
            PERSISTENT.Delta_old = DELTA_OLD;
            return
        end
        
        if newnrm < epsilon
            SOL = Y_guess-dY;
            break;
        else
            Y_guess = Y_guess-dY;
            gY = g(t,Y_guess);
            ODE_CALLS = ODE_CALLS + 1;
            oldnrm = newnrm;
        end
    end
    if iter == itmax
        %disp(['Maximum no of iterations reached, hh, h = ', num2str(dtCellVec(end))])
        SOL = NaN(size(Y_TYPICAL));
        DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
        %PERSISTENT=[{FAC} {DELTA_OLD}];
        PERSISTENT.Fac = FAC;
        PERSISTENT.Delta_old = DELTA_OLD;
        return
    end
%     if t == h
%         DELTA_OLD = [zeros(size(Y_old)), Y_old - SOL, zeros(size(Y_old))];
%     else
%         DELTA_OLD(:,end) = Y_old - SOL;
%     end
    DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), Y_old - SOL];
    %FAC=[];
else  % step > 2
    % Setting up parameters for the BDF2 iteration
    gamma = DT(end)/DT(end-1);
    beta = (gamma + 1)/(2*gamma + 1);
    alpha2 = -gamma^2/(2*gamma+1);
    alpha1 = 1-alpha2;
    
    % setting up parameters for the second order guess
    delta = 1+DT(end-2)/DT(end-1);    % 1 + h(n-2)/h(n-1)
    a2 = gamma*(gamma+delta)/(1-delta);
    a3 = gamma*(gamma+1)/(delta*(delta-1));
    a1 = 1-a2-a3;
    
    Y_old = Y(:,end);
    Y_ancient = Y(:,end-1);
    Y0 = a1*Y_old+a2*Y_ancient + a3*Y(:,end-2);
    Y0_dot = (1/(beta*DT(end)))*(Y0 - alpha1*Y_old - alpha2*Y_ancient);
    
    Delta_guess = DELTA_OLD(:,end) + gamma*(DELTA_OLD(:,end)-DELTA_OLD(:,end-1));
    
    %        g = @(t,Delta) bsxfun(@minus,Delta, ...
    
    %            beta*dtCellVec(Nhh*(sysIndex-1)+hh_index)* ...
    
    %            bsxfun(@minus,Y0_dot,ode_cell(t, bsxfun(@minus,Y0,Delta),frac,cai)));
    
    f = @(t,SOL)feval(ODE_FUN,t, SOL, params{:});
    
    g = @(t,Delta) (Delta - beta*DT(end)*(Y0_dot - feval(ODE_FUN,t, Y0-Delta, params{:})));
    
    %        g = @(t, Y_new) bsxfun(@minus,bsxfun(@minus,Y_new,alpha1*Y_old),bsxfun(@plus,alpha2*Y_ancient,...
    
    %            beta*dtCellVec(Nhh*(sysIndex-1)+hh_index)*ode_cell(t,Y_new, frac, cai)));
    
    
    
    % Prepare for the modified Newton iteration

    gY = g(t,Delta_guess);
    fY = f(t, Y0-Delta_guess);
    ODE_CALLS = ODE_CALLS+1;
    
    [J, FAC,~,nf] = numjac(f,t,Y0-Delta_guess,fY,{1e-6,Y_TYPICAL}, FAC,0);
    NJ_CALLS = NJ_CALLS+1;
    NJ_ODE_CALLS = NJ_ODE_CALLS + nf;
        
    [LJ, UJ] = lu(eye(size(J))-beta*DT(end)*J);
    
    % Solve via iteration
    while (iter < itmax)
        iter  = iter + 1;
        
        dDelta = UJ\(LJ\gY);
                
        newnrm = norm(dDelta,Inf);
        if iter > 1 && newnrm >= 0.9*oldnrm
            disp(['Slowly convergent Newton method, h = ', num2str(DT(end))])
            SOL = NaN(size(Y_TYPICAL));
            DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
            %PERSISTENT=[{FAC} {DELTA_OLD}];
            PERSISTENT.Fac = FAC;
            PERSISTENT.Delta_old = DELTA_OLD;
            return
        end
        
        if newnrm < epsilon
            SOL = Y0 - (Delta_guess - dDelta);
            %yCellVec(:,Nhh*(sysIndex-1)+hh_index + 1) = Y_guess-dY;
            break;
        else
            Delta_guess = Delta_guess-dDelta;
            gY = g(t,Delta_guess);
            ODE_CALLS = ODE_CALLS + 1;
            oldnrm = newnrm;
        end
    end
    
    if (iter == itmax)
        %disp(['Maximum no of iterations reached, hh, h = ', num2str(dtCellVec(end))])
        SOL = NaN(size(Y_TYPICAL));
        DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
        %PERSISTENT=[{FAC} {DELTA_OLD}];
        PERSISTENT.Fac = FAC;
        PERSISTENT.Delta_old = DELTA_OLD;
        return
    end
    DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), Y0-SOL];
end

%PERSISTENT=[{FAC} {DELTA_OLD}];  
PERSISTENT.Fac = FAC;
PERSISTENT.Delta_old = DELTA_OLD;

end
