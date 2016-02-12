% MAX_EEST :
%   the maximum error value estimated by other controllers,
%   valid in multirate
%
function [H_MACRO, H_MICRO STEP_REJECTED, PERSISTENT] = ec_h211b( DT, ...
    SYS_INDEX,  E_EST,  PERSISTENT)

% The limit for how much the steps are allowed to grow

rho = 0.90;     % Safety factor for tolerance limit

dtMax = Inf;    % Maximum step size

TOL = 1;        % Tolerance level for error

eEstVec = PERSISTENT.eEstVec; % Estimations of last errors

rhofac = PERSISTENT.rhofac; % Söderlind's limiter

STEP_REJECTED = false;

deuflhard = false;

q = 2;

% both constants used for micro time step prediction
hmax = 1e-2;
mmax = 2;


%     H_MACRO = DT(end)/2;
%     H_MICRO = PERSISTENT.H_MICRO/2;
%
%     %     if (any(isnan(in(1:erkSize))))
%     %         nanItersErk = nanItersErk+1;
%     %     end
%     %
%     STEP_REJECTED = true;



if SYS_INDEX > 1  % for later steps, find the error estimate
    
    if(SYS_INDEX == 2)
        
        rhofac = DT(end)/DT(end-1);
        
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
        
        
        
        %         Deuflhard's controller
        if (deuflhard)
            
            dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
                (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
                (rho*TOL/eEstVec(1))^(c/p);     % Differenting part
            
            
        else %         H211b controller by Söderlind
            
            rhofac = (rho*TOL/eEstVec(3))^(0.25/p) * ...
                ((rho*TOL)/eEstVec(2))^(0.25/p)* ...
                (rhofac)^(-0.25);
            
            kappa = 1;
            
            rhofac = kappa*atan((rhofac-1)/kappa)+1;
            
            dtSuggest = DT(end)*rhofac;
        end
        
        H_MACRO = min([q*DT(end), dtMax, dtSuggest]);
        H_MICRO = min(hmax, PERSISTENT.h / max(0.2, 1.25*eEstVec(3)^(-p)));
        PERSISTENT.M = max(fix(PERSISTENT.M/1.5),max(mmax,fix(H_MACRO/H_MICRO)));
        
        
        STEP_REJECTED = false;
        
    else
        
        %         Classical controller
        if (deuflhard)
            dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
                (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
                (rho*TOL/eEstVec(1))^(c/p);
            %                 %         H211b controller by Söderlind
            %                 dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(0.25/p) * ...
            %                                ((rho*TOL)/eEstVec(2))^(0.25/p)* ...
            %                                (dt(end)/dt(end-1))^(-0.25);
        else
            dtSuggest = DT(end)*(rho*TOL/eEstVec(3))^(1/p);
        end
        % rhofac = dt(end)/dt(end-1);
        
        
        
        %fold = (rho*TOL/eEstVec(3))^(1/p);
        
        H_MACRO = min([q*DT(end), dtMax, dtSuggest]);
        H_MICRO = PERSISTENT.h*max(.1,.8*eEstVec(3)^p);
       
        PERSISTENT.M=min(fix(1.15*PERSISTENT.M),...
            max(mmax,fix(H_MACRO/H_MICRO)));
        
        eEstVec = circshift(eEstVec',1)';
        
        STEP_REJECTED = true;
        
    end
    
    
else % For the first step, use Euler formula with no error estimation
    
    H_MACRO = DT(end);
    PERSISTENT.M = mmax;
    eEstVec = [NaN rho*TOL rho*TOL];
    
end
PERSISTENT.eEstVec = eEstVec;
PERSISTENT.rhofac = rhofac;

H_MICRO=H_MACRO/PERSISTENT.M;
PERSISTENT.h = H_MICRO;

end