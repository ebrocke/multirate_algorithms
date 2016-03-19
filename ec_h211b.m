
function [H STEP_REJECTED PERSISTENT] = ec_h211b( DT, ...
    E_EST,  PERSISTENT)

% The limit for how much the steps are allowed to grow

rho = 0.90;     % Safety factor for tolerance limit

dtMax = Inf;    % Maximum step size

TOL = 1;        % Tolerance level for error


STEP_REJECTED = false;

%deuflhard = false;

q = 2;
if isinf(E_EST)
    H = DT(end)/2;
    STEP_REJECTED = true;
    return
end

if PERSISTENT.init > 0  % for later steps, find the error estimate
    eEstVec = PERSISTENT.eEstVec; % Estimations of last errors
    
    if(PERSISTENT.init == 1)
        rhofac = DT(end)/DT(end-1);
    else
        rhofac = PERSISTENT.rhofac; % Söderlind's limiter

    end
    
    
    
    
    % Once we have an error estimation, define the optimal timestep
    
    % for the current time. Order 2 method => 3
    
    
    
    eEstVec = circshift(eEstVec',2)';
    
    eEstVec(3) = E_EST;
    
    
    
    % Using notation from Hairer, Deuflhard uses a different
    
    % notation, but does not specify the constants...
    
    p = 2; % Order of the method+1
    
    
    a = 0.7;    % 1
    b = 0.4;    % 0
    c = 0;
    
    %             if step_rejected
    
    %                 dtSuggest = dt(end)*(rho*TOL/eEstVec(3))^(1/p);
    
    %             else
    
    %                 dtSuggest = dt(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
    
    %                   (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
    
    %                   (rho*TOL/eEstVec(1))^(c/p);     % Differenting part
    
    %             end
    
    
    
    %             diffTvalues = [q*dt(end), dtMax, dtSuggest];
    
    %             dt_next = min(diffTvalues);
    
    if E_EST < 1
        
        % if the error is smaller than the tolerance, step is
        
        % accepted and the optimal time step for the current step
        
        % will be used as the next time step.
        
        
        
        %         %         Deuflhard's controller
        %         if (deuflhard)
        %
        %             dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
        %                 (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
        %                 (rho*TOL/eEstVec(1))^(c/p);     % Differenting part
        %
        %
        %         else %         H211b controller by Söderlind
        
        rhofac = (rho*TOL/eEstVec(3))^(0.25/p) * ...
            ((rho*TOL)/eEstVec(2))^(0.25/p)* ...
            (rhofac)^(-0.25);
        
        kappa = 1;
        
        rhofac = kappa*atan((rhofac-1)/kappa)+1;
        
        dtSuggest = DT(end)*rhofac;
        %        end
        
        H = min([q*DT(end), dtMax, dtSuggest]);
        
        PERSISTENT.init = 2;
        STEP_REJECTED = false;
        
    else
        
        %         %         Classical controller
        %         if (deuflhard)
        %             dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
        %                 (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
        %                 (rho*TOL/eEstVec(1))^(c/p);
        %             %                 %         H211b controller by Söderlind
        %             %                 dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(0.25/p) * ...
        %             %                                ((rho*TOL)/eEstVec(2))^(0.25/p)* ...
        %             %                                (dt(end)/dt(end-1))^(-0.25);
        %         else
        dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(1/p);
        %        end
        % rhofac = dt(end)/dt(end-1);
        
        %fold = (rho*TOL/eEstVec(3))^(1/p);
        
        H = min([q*DT(end), dtMax, dtSuggest]);
        
        eEstVec = circshift(eEstVec',1)';
        
        STEP_REJECTED = true;
        
    end
    
    
else % For the first step, use Euler formula with no error estimation
    PERSISTENT.init = 1;
    eEstVec = [NaN rho*TOL rho*TOL];
    rhofac = 0;
    H = 1e-5;
    
end
PERSISTENT.h = H;
PERSISTENT.eEstVec = eEstVec;
PERSISTENT.rhofac = rhofac;

end