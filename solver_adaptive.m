function [t_erk, y_erk, t_cell, y_cell, stats_erk, stats_cell]=...
    solver_adaptive(t, ...
    initVals, ...
    relTol, ...
    erkSize, ...
    erkSys, cellSys,...
    multirate)

global MODE 


tEnd = t(2);     % Ending time
t_ = t(1);         % Current time


step_rejected_ = false;

yTypical = importdata('yTypicalSolution.txt',';'); % Typical solution

% memory allocation constants
% we enlarge by k1 if full
% not efficient for multirate
k1 = 500000;
baseTol=1e-4;
mem_size=k1*2^(log10(relTol/baseTol));

%% Persistent variables can be considered as private working variables
% statistics
P_ERK.stats.acceptedIter = 0;
P_ERK.stats.refinedIter = 0;
P_ERK.stats.rejectedIter = 0;
P_ERK.stats.numjac = 0;
P_ERK.stats.n_ode_numjac = 0;
P_ERK.stats.n_ode_iter = 0;
% solver working variables
P_ERK.solver.init = true;
P_ERK.solver.yTypical = yTypical(1:erkSize);
% controller working variables
P_ERK.controller.eEst = 0;
P_ERK.controller.init = 0;
P_ERK.controller.h = 1e-5;
% system working variables
P_ERK.sys.t_exch = zeros(1,MODE); %store last exchanged values
P_ERK.sys.y_exch = zeros(MODE,2); % in a reversed order of time
% handles that provides flexibility in solving the system
P_ERK.sys.method_hdl = erkSys{1}; % BDF_DEF BDF_AS RK4 CrankNicStagg; 
P_ERK.sys.exch_hdl = erkSys{2}; % system function that provides exchanged variables
P_ERK.sys.isolver_hdl = erkSys{3};% how to solve interval (used in multirate)
P_ERK.sys.ode_hdl = @ode_erk; % handle to the system ode functions
% store system solution
P_ERK.sol.y = zeros(erkSize,mem_size);
P_ERK.sol.dt = zeros(1,mem_size);
P_ERK.sol.y(:,1) = initVals(1:erkSize)';


P_CELL.solver.init = true;
P_CELL.solver.yTypical = yTypical(erkSize+1:end);
P_CELL.stats.acceptedIter = 0;
P_CELL.stats.refinedIter = 0;
P_CELL.stats.rejectedIter = 0;
P_CELL.stats.numjac = 0;
P_CELL.stats.n_ode_numjac = 0;
P_CELL.stats.n_ode_iter = 0;
P_CELL.sol.y = zeros(length(initVals)-erkSize,mem_size);
P_CELL.sol.dt = zeros(1,mem_size);
P_CELL.sol.y(:,1) = initVals(erkSize+1:end)';
P_CELL.controller.init = 0;
P_CELL.controller.eEst = 0;
P_CELL.controller.h = 1e-5;
P_CELL.sys.method_hdl = cellSys{1};
P_CELL.sys.exch_hdl = cellSys{2};
P_CELL.sys.isolver_hdl = cellSys{3};
P_CELL.sys.ode_hdl = @ode_cell;
P_CELL.sys.t_exch = zeros(1,MODE);
P_CELL.sys.y_exch = zeros(MODE,2);
P_CELL.sys.h = 1e-5;
P_CELL.sys.m = 1;
P_CELL.sys.eEst  = 0;

PERS.ERK = P_ERK;
PERS.CELL = P_CELL;
PERS.yTypical = yTypical;

yTypical = abs(yTypical);
H_ = 1e-5;

% We save last three system solutions
% is the error estimation
y_ = zeros(length(initVals),3);
dt_ = zeros(1,3);
y_(:,end) = initVals;
dt_(end) = H_;


if multirate
   ysize = erkSize;
   yTypical = P_ERK.solver.yTypical;
else
    ysize = length(yTypical);
end
h_ = 1e-5;
H_ = 1e-5;

while t_ < tEnd
    
    ii_ = PERS.ERK.stats.acceptedIter;
    if(rem(ii_, 1000)==0) % for displaying progress
        
         toc, t_, H_, h_
        
    end

%     if (h_ < 5e-5)
%         stop = 0 ;
%         [H_ h_ step_rejected_]
%     end

    
    % Calculate the solution
    [sol_,  PERS] = solve_sys(...
        [t_ t_ + H_], relTol, step_rejected_, PERS);
    %dbstop if warning
    % calculate the error
    [eEst_, eI] = ee_skelboe2000(...
        [y_(1:ysize,:) sol_(1:ysize)], dt_, relTol, yTypical);
    PERS.CELL.sys.eEst = eEst_;
    % calculate the macro H_ and micro h_ time steps
     % we do not change macro time step if micro time step was not fine
     % enough
%     if any(isnan(sol_(erkSize+1:end)))
%         H_ = dt_(end);
%         step_rejected_ = true;
%     else
%    [H_,  step_rejected_, P_ERK.controller] = ec_h211b(dt_, ...
%        max(PERS.CELL.controller.eEst, eEst_),...
%        P_ERK.controller);
%    end
   % k =  1/(H_/PERS.CELL.controller.eEst)^2;
    %[eEst_ eEst_*k]
    if ~isnan(eEst_)
        eEst_ = max(eEst_,PERS.CELL.controller.eEst);
    end
    [H_,  step_rejected_, P_ERK.controller] = ec_h211b(dt_, ...
        eEst_,...
        P_ERK.controller);

    
%     if (multirate) %ec_h211b_hmicro
%         [h_ P_CELL.controller] = ec_classical(H_,...
%             max(PERS.CELL.controller.eEst, eEst_),...
%             P_CELL.controller);
%         PERS.CELL.sys.h = h_; % suggested micro time step
%         PERS.CELL.sys.m = P_CELL.controller.m; % suggensted n of intervals
% 
%     %    [H_ h_ step_rejected_]
%     end
%        [H_ h_ step_rejected_ P_CELL.controller] = ec_classical_comb(...
%            dt_(end),...
%             max(PERS.CELL.controller.eEst, eEst_),...
%             P_CELL.controller);
%         PERS.CELL.sys.h = h_;
%         PERS.CELL.sys.m = P_CELL.controller.m;

    if (~step_rejected_)
        t_ = t_ + dt_(end);
        dt_=circshift(dt_,[0,-1]);
        y_=circshift(y_,[0,-1]);
        y_(:,end) = sol_;
    end
    dt_(end) = H_;
    % Ensure that we dont leave the interval
    if t_ + H_ > tEnd  
        H_ = tEnd-t_;
    end
    
    % add check up for memory allocation in the solution vectors
    
end

% fill in return values
i_ = PERS.ERK.stats.acceptedIter;
y_erk = PERS.ERK.sol.y(:, 1:i_+1);
t_erk = [0 cumsum(PERS.ERK.sol.dt(1:i_))];
stats_erk = PERS.ERK.stats;

i_ = PERS.CELL.stats.acceptedIter;
y_cell = PERS.CELL.sol.y(:, 1:i_+1);
t_cell = [0 cumsum(PERS.CELL.sol.dt(1:i_))];
stats_cell = PERS.CELL.stats;

% 
 fprintf('%-20s \t %-8s \t %-8s\n',sprintf('System'), 'Erk','HH')
% 
 fprintf('%-20s \t %-8d \t %-8d\n','numjac calls', ...
     PERS.ERK.stats.numjac, PERS.CELL.stats.numjac)
% 
 fprintf('%-20s \t %-8d \t %-8d\n','n ode by numjac', ...
     PERS.ERK.stats.n_ode_numjac, PERS.CELL.stats.n_ode_numjac)
% 
 fprintf('%-20s \t %-8d \t %-8d\n','n ode by iterations',...
     PERS.ERK.stats.n_ode_iter, PERS.CELL.stats.n_ode_iter)
%
  fprintf('%-20s \t %-8d \t %-8d\n','accepted calls', ...
     PERS.ERK.stats.acceptedIter, PERS.CELL.stats.acceptedIter)
%
  fprintf('%-20s \t %-8d \t %-8d\n','refined calls', ...
     PERS.ERK.stats.refinedIter, PERS.CELL.stats.refinedIter)
 %
  fprintf('%-20s \t %-8d \t %-8d\n','rejected calls', ...
     PERS.ERK.stats.rejectedIter, PERS.CELL.stats.rejectedIter)
end
