function [SOL, NJ_CALLS, NJ_ODE_CALLS, ODE_CALLS, PERSISTENT] = BDF2_DEF(T,...
    Y, DT, ODE_FUN, ODE_PARAMS, relTol, PERSISTENT)

Y_TYPICAL = PERSISTENT.yTypical;
STEP_REJECTED = PERSISTENT.step_rejected;
if PERSISTENT.init
    FAC=[];
    DELTA_OLD = zeros(size(Y_TYPICAL,1), 3);
    PERSISTENT.init = false;
else
    FAC = PERSISTENT.Fac;
    DELTA_OLD = PERSISTENT.Delta_old;
    NEW_JACOBIAN = PERSISTENT.New_Jac;
    J = PERSISTENT.J;
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

epsilon = 0.01*relTol*norm(Y_TYPICAL,Inf);
itmax = 5;
iter = 0;
h = T(2)-T(1);
t = T(2); % one step


if size(Y,2) <=2 %t == h || size(Y,2) == 2  % step 1 and 2
    
    % Implicit Euler step
    Y_old = Y(:,end);    
    g = @(t,Y_new) Y_new-Y_old-DT(end)*feval(ODE_FUN,t,Y_new, params{:});
    
    %        g = @(t,Y_new) bsxfun(@minus,bsxfun(@minus,Y_new,Y_old),...
    
    %            dtCellVec(end)*ode_cell(t,Y_new, frac, cai));
    
    % Constant interpolation
    Y_guess = Y_old;
    
     % Prepare for the modified newton iteration
    gY = g(t,Y_guess);
    ODE_CALLS = ODE_CALLS + 1;
    
    [J, FAC, ~, nf] = numjac(g,t,Y_guess,gY,{1e-6,Y_TYPICAL}, FAC,0);
    NJ_CALLS = NJ_CALLS + 1;
    NJ_ODE_CALLS = NJ_ODE_CALLS + nf;
    [LJ, UJ] = lu(J);
 
    % Solve via iteration
    while (iter < itmax)
        
        iter = iter+1;
        dY = UJ\(LJ\gY);
        
        newnrm = norm(dY,Inf);
        if iter > 1 && newnrm >= 0.9*oldnrm
            %disp(['Slowly convergent Newton method, hh, h = ', num2str(dtCellVec(end))])
            SOL = NaN(size(Y_TYPICAL));
            DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
            %PERSISTENT=[{FAC} {DELTA_OLD} {NEW_JACOBIAN} {J}];
            PERSISTENT.Fac = FAC;
            PERSISTENT.Delta_old = DELTA_OLD;
            PERSISTENT.New_Jac = NEW_JACOBIAN;
            PERSISTENT.J = J;
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
    
    if (iter == itmax)
        %disp(['Maximum no of iterations reached, hh, h = ', num2str(dtCellVec(end))])
        SOL = NaN(size(Y_TYPICAL));
        DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
        %PERSISTENT=[{FAC} {DELTA_OLD} {NEW_JACOBIAN} {J}];
        PERSISTENT.Fac = FAC;
        PERSISTENT.Delta_old = DELTA_OLD;
        PERSISTENT.New_Jac = NEW_JACOBIAN;
        PERSISTENT.J = J;
        return
    end
    %if t == h
    %    DELTA_OLD = [zeros(size(Y_old)), Y_old - SOL, zeros(size(Y_old))];
    %else
    DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), Y_old - SOL];
    %end
    
    NEW_JACOBIAN = true;
    %njac = 50;
    %FAC = [];
    
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

    Delta_guess = DELTA_OLD(:,end) + ...
        gamma*(DELTA_OLD(:,end)-DELTA_OLD(:,end-1));

%        g = @(t,Delta) bsxfun(@minus,Delta, ...

%            beta*dtCellVec(Nhh*(sysIndex-1)+hh_index)* ...

%            bsxfun(@minus,Y0_dot,ode_cell(t, bsxfun(@minus,Y0,Delta),frac,cai)));

    f = @(t,SOL)feval(ODE_FUN, t, SOL, params{:});

    g = @(t,Delta) (Delta - beta*DT(end)*...
        (Y0_dot - feval(ODE_FUN, t, Y0-Delta, params{:})));

%        g = @(t, Y_new) bsxfun(@minus,bsxfun(@minus,Y_new,alpha1*Y_old),bsxfun(@plus,alpha2*Y_ancient,...

%            beta*dtCellVec(Nhh*(sysIndex-1)+hh_index)*ode_cell(t,Y_new, frac, cai)));



    % Prepare for the modified newton iteration
    state = true;
    gY_old = g(t,Delta_guess);
    fY = f(t, Y0-Delta_guess);
    ODE_CALLS = ODE_CALLS+1;

    while state
        if NEW_JACOBIAN
            [J, FAC, ~, nf] = numjac(f,t,Y0-Delta_guess,fY,{1e-6,Y_TYPICAL}, FAC,0);
            NJ_CALLS = NJ_CALLS+1;
            NJ_ODE_CALLS = NJ_ODE_CALLS+nf;
            %disp(['Number of jac-evals: ',num2str(njac)])
        end
        
        [LJ, UJ] = lu(eye(size(J))-beta*DT(end)*J);
        gY = gY_old;
        Delta_g = Delta_guess;
        
        iter = 0;
        % Solve via iteration
        while (iter < itmax)

            iter = iter + 1;
            dDelta = UJ\(LJ\gY);

            newnrm = norm(dDelta,Inf);
            if iter > 1 && newnrm >= 0.9*oldnrm
                if NEW_JACOBIAN
                    %disp(['Slowly convergent Newton method, hh, h = ', num2str(dtCellVec(end))])
                    SOL = NaN(size(Y_TYPICAL));
                    DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
                    %PERSISTENT=[{FAC} {DELTA_OLD} {NEW_JACOBIAN} {J}];
                    PERSISTENT.Fac = FAC;
                    PERSISTENT.Delta_old = DELTA_OLD;
                    PERSISTENT.New_Jac = NEW_JACOBIAN;
                    PERSISTENT.J = J;
                    return
                else
                    NEW_JACOBIAN = true;
                    break
                end
            end
            if newnrm < epsilon
                SOL = Y0 - (Delta_g - dDelta);
                %yCellVec(:,Nhh*(sysIndex-1)+hh_index + 1) = Y_guess-dY;
                state = false;
                break;
            else
                Delta_g = Delta_g-dDelta;
                gY = g(t,Delta_g);
                ODE_CALLS = ODE_CALLS + 1;
                oldnrm = newnrm;
            end
        end
        
        if (iter == itmax)
            if NEW_JACOBIAN
                %disp(['Maximum no of iterations reached, hh, h = ', num2str(dtCellVec(end))])
                SOL = NaN(size(Y_TYPICAL));
                DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), SOL];
                %PERSISTENT=[{FAC} {DELTA_OLD} {NEW_JACOBIAN} {J}];
                PERSISTENT.Fac = FAC;
                PERSISTENT.Delta_old = DELTA_OLD;
                PERSISTENT.New_Jac = NEW_JACOBIAN;
                PERSISTENT.J = J;
                return
            else
                NEW_JACOBIAN = true;
            end
        end
    end
    DELTA_OLD = [DELTA_OLD(:,end-1), DELTA_OLD(:,end), Y0-SOL];
    NEW_JACOBIAN = false;
end

%PERSISTENT=[{FAC} {DELTA_OLD} {NEW_JACOBIAN} {J}];
PERSISTENT.Fac = FAC;
PERSISTENT.Delta_old = DELTA_OLD;
PERSISTENT.New_Jac = NEW_JACOBIAN;
PERSISTENT.J = J;

end