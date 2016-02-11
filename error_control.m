function [NEW_DT, STEP_REJECTED, PERSISTENT] = error_control(Y, DT, ...
    SYS_INDEX, RELTOL, YTYPICAL, PERSISTENT, err)
% TODO: the external controller should be considered


% The limit for how much the steps are allowed to grow

rho = 0.90;     % Safety factor for tolerance limit

dtMax = Inf;    % Maximum step size

TOL = 1;        % Tolerance level for error

eEstVec = PERSISTENT.eEstVec; % Estimations of last errors

rhofac = PERSISTENT.rhofac; % Söderlind's limiter

STEP_REJECTED = false;

eEst = NaN;

deuflhard = false;

q = 2;

if any(isnan(Y (:,end))) % the last element is the solution to evaluate
    
    NEW_DT = DT(end)/2;
    
    %     if (any(isnan(in(1:erkSize))))
    %         nanItersErk = nanItersErk+1;
    %     end
    %
    STEP_REJECTED = true;
   
    
else
    
    if SYS_INDEX > 1  % for later steps, find the error estimate
        
        if(SYS_INDEX == 2)
            
            rhofac = DT(end)/DT(end-1);
        end
        
        [eEst, ~] = error_estimate(Y, DT, SYS_INDEX, RELTOL, YTYPICAL);
        
        
        % Once we have an error estimation, define the optimal timestep
        
        % for the current time. Order 2 method => 3
        
        
        
        eEstVec = circshift(eEstVec',2)';

        eEstVec(3) = max(eEst,err);
        
        
        
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
        
        if eEst < 1
            
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
            
            NEW_DT = min([q*DT(end), dtMax, dtSuggest]);
            
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
            
            NEW_DT = min([q*DT(end), dtMax, dtSuggest]);
            
            eEstVec = circshift(eEstVec',1)';
            
            STEP_REJECTED = true;
            
        end
        
        
    else % For the first step, use Euler formula with no error estimation
        
        NEW_DT = DT(end);
        eEstVec = [NaN rho*TOL rho*TOL];
        
    end
    PERSISTENT.eEstVec = eEstVec;
    PERSISTENT.rhofac = rhofac;
    PERSISTENT.err = eEst;
end
end