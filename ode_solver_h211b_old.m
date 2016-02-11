function [tInt, y]=ode_solver_h211b(solveSys,tInt, ...
    initVals, relTol, erkSize, erkSys, cellSys)

%%%% initializing global variables
%%%% Function counters and controll parameters

global nanIters nanItersErk refineIters acceptedErrorIndex rejectedErrorIndex
global dt_ratio y_ratio RATIO
global dt sysIndex stats
% Inital time step (used by adaptive solver)
dt = 1e-5;
RATIO = 1;
%global a b

%%%% Defining time intervals and step sizes

tStart = tInt(1);   % Starting time

tEnd = tInt(2);     % Ending time

t = tStart;         % Current time



%Nsys = 1; 

%hhIter = zeros(1,Nsys);

%erkIter = zeros(1,Nsys);

%yCellVec = zeros(length(initVals(erkSize+1:end)),Nsys);

%yCellVec(:,1) = initVals(erkSize+1:end);

%yErkVec = zeros(length(initVals(1:erkSize)),Nsys);

%yErkVec(:,1) = initVals(1:erkSize);



step_rejected = false;



%%%% Defining parameters for the adaptive solver

q = 2;          % The limit for how much the steps are allowed to grow

rho = 0.90;     % Safety factor for tolerance limit

dtMax = Inf;    % Maximum step size

TOL = 1;        % Tolerance level for error



yTypical = importdata('yTypicalSolution.txt',';'); % Typical solution

eEstVec = [NaN rho*TOL rho*TOL]; % Estimations of last errors



% Persistent now handled in a more concise way

P_ERK.init = true;

P_ERK.yTypical = yTypical(1:erkSize);

P_CELL.init = true;

P_CELL.yTypical = yTypical(erkSize+1:end);

PERS.ERK = P_ERK;

PERS.CELL = P_CELL;



% Set up initial values etc.

y = zeros(length(initVals),1);

y(:,1) = initVals;

yTypical = abs(yTypical);
Ypn = zeros(length(initVals),1);
Yp2n = zeros(length(initVals),1);





while t < tEnd

    sysIndex = sysIndex + 1;

    if(rem(sysIndex,1000)==0) % for displaying progress

        sysIndex, toc, t

    end

   
    % Calculate the solution

    [out, stat, PERS]=solveSys([t t + dt(end)], y, erkSize,...
        erkSys, cellSys, relTol, PERS, step_rejected);

    stats = stats+stat;

    % Update time and timestep

    if any(isnan(out))

        dt(end) = dt(end)/2;

        sysIndex = sysIndex - 1;

        nanIters = nanIters + 1;
        if (any(isnan(out(1:erkSize))))
            nanItersErk = nanItersErk+1;
        end

        step_rejected = true;

    else

        if sysIndex == 1  % For the first step, use Euler formula with no error estimation

            t = t + dt(end);

            y(:,sysIndex+1) = out;

            dt(end+1) = dt(end);

        else % for later steps, find the error estimate

            if sysIndex == 2 % For the second step, use BDF2 with mode 2

                gamma = dt(end)/dt(end-1);  % = h(n)/h(n-1)
                Ypn(1:erkSize) = y(1:erkSize,sysIndex)...
                    +gamma*(y(1:erkSize,sysIndex)...
                    -y(1:erkSize,sysIndex-1));
                
                gamma = dt_ratio(end)/dt_ratio(end-1);
                Ypn(erkSize+1:end) = y_ratio(:,end-1)...
                    +gamma*(y_ratio(:,end-1)...
                    -y_ratio(:,end-2));

                

                % Should it be ABS (Y_typ) or only (Y_typ)?????

                [eEst, eIndex] = max(abs(out-Ypn)./(relTol*(abs(out)+yTypical)));

                % Söderlind's limiter

                rhofac = dt(end)/dt(end-1);

                

           else % for later steps, use BDF2 with mode 3

               gamma = dt(end)/dt(end-1);
               delta = 1+dt(end-2)/dt(end-1);    % 1 + h(n-2)/h(n-1)
               a2 = gamma*(gamma+delta)/(1-delta);
               a3 = gamma*(gamma+1)/(delta*(delta-1));
               a1 = 1-a2-a3;
               Yp2n(1:erkSize) = a1*y(1:erkSize,sysIndex)...
                   +a2*y(1:erkSize,sysIndex-1)...
                   +a3*y(1:erkSize,sysIndex-2);

               
               gamma = dt_ratio(end)/dt_ratio(end-1);
               delta = 1+dt_ratio(end-2)/dt_ratio(end-1);    % 1 + h(n-2)/h(n-1)
               a2 = gamma*(gamma+delta)/(1-delta);
               a3 = gamma*(gamma+1)/(delta*(delta-1));
               a1 = 1-a2-a3;
               Yp2n(erkSize+1:end) = a1*y_ratio(:,end-1)...
                   +a2*y_ratio(:,end-2)...
                   +a3*y_ratio(:,end-3);

               %This has an explanation...

               [eEst, eIndex] = max(abs(out-Yp2n)./(relTol*(abs(out)+yTypical)));

           end

            

            % Once we have an error estimation, define the optimal timestep

            % for the current time. Order 2 method => 3

            

            eEstVec = circshift(eEstVec',2)';

            eEstVec(3) = eEst;

            

            % Using notation from Hairer, Deuflhard uses a different

            % notation, but does not specify the constants...

            p = 2; % Order of the method+1

            a = 1;
            b = 0;
            %a = 0.7; 
            
            %b = 0.4; 

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

                t = t + dt(end);

%         Deuflhard's controller

%               dtSuggest = dt(end)*(rho*TOL/eEstVec(3))^(a/p) * ... %Proportional part
%                   (eEstVec(2)/(rho*TOL))^(b/p)* ... % Integrating part
%                   (rho*TOL/eEstVec(1))^(c/p);     % Differenting part
%         H211b controller by Söderlind

                rhofac = (rho*TOL/eEstVec(3))^(0.25/p) * ... 
                  ((rho*TOL)/eEstVec(2))^(0.25/p)* ... 
                  (rhofac)^(-0.25);     

                kappa = 1;

                rhofac = kappa*atan((rhofac-1)/kappa)+1;

                dtSuggest = dt(end)*rhofac;

                dt(end+1) = min([q*dt(end), dtMax, dtSuggest]);

                y(:,sysIndex+1) = out;

                acceptedErrorIndex(end+1) = eIndex;

                step_rejected = false;

            else

%         Classical controller

                 dtSuggest = dt(end)*(rho*TOL/eEstVec(3))^(1/p);

                 rhofac = dt(end)/dt(end-1);

%         H211b controller by Söderlind

%                 dtSuggest = dt(end)*(rho*TOL/eEstVec(3))^(0.25/p) * ... 
%                   ((rho*TOL)/eEstVec(2))^(0.25/p)* ... 
%                   (dt(end)/dt(end-1))^(-0.25);     

                fold = (rho*TOL/eEstVec(3))^(1/p);

                dt(end) = min([q*dt(end), dtMax, dtSuggest]);

                refineIters = refineIters + 1;

                sysIndex = sysIndex-1;

                eEstVec = circshift(eEstVec',1)';

                rejectedErrorIndex(end+1) = eIndex;

                step_rejected = true;

            end

        end

        

        if t + dt(end) > tEnd  % Ensure that we dont leave the interval

            dt(end) = tEnd-t;

        end

    end

end

%if(~sum(erkIter)), erkIter = length(dt); end

%if(~sum(hhIter)), hhIter = Nsys*Nhh; end



tInt = [0 cumsum(dt)];

%y = y(:,1:length(dt)+1);

fprintf('%d accepted steps\n',sysIndex)

fprintf('%d steps had to be refined and taken again\n',refineIters)

fprintf('%d steps (%d erk component) resulted in NaN and had to be taken again\n',nanIters, nanItersErk)

fprintf('%-20s \t %-8s \t %-8s\n',sprintf('System'), 'Erk','HH')

% fprintf('%-20s \t %-8d \t %-8d\n','Solver calls',sum(erkIter),sum(hhIter))

% fprintf('%-20s \t %-8d \t %-8d\n','Failed solver calls',...

%     sum(erkIter)-length(erkIter),sum(hhIter)-length(hhIter))

%fprintf('%-20s \t %-8d \t %-8d\n','ode calls',ode_erk_iter,ode_cell_iter)

%fprintf('%-20s \t %-8d \t %-8d \t %-8d\n ','numjac calls',numjac_calls_erk - ode_erk_iter,...

%    numjac_calls_cell, test_numjac)

fprintf('%-20s \t %-8d \t %-8d\n','numjac calls', stats(1,1), stats(2,1))

fprintf('%-20s \t %-8d \t %-8d\n','ode calls by numjac', stats(1,2), stats(2,2))

fprintf('%-20s \t %-8d \t %-8d\n','ode calls by iterations', stats(1,3), stats(2,3))



end