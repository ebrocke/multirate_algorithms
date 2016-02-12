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
k1 = 50000;
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
P_ERK.controller.eEstVec = 0;
P_ERK.controller.rhofac = 0;
P_ERK.controller.eEst = 0;
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
P_CELL.controller.eEst = 0;
P_CELL.controller.eEstVec = 0;
P_CELL.controller.rhofac = 0;
P_CELL.controller.optimal_dt = 1e-5;
P_CELL.sys.method_hdl = cellSys{1};
P_CELL.sys.exch_hdl = cellSys{2};
P_CELL.sys.isolver_hdl = cellSys{3};
P_CELL.sys.ode_hdl = @ode_cell;
P_CELL.sys.t_exch = zeros(1,MODE);
P_CELL.sys.y_exch = zeros(MODE,2);


PERS.ERK = P_ERK;
PERS.CELL = P_CELL;
PERS.yTypical = yTypical;

yTypical = abs(yTypical);
dt_optimal_ = 1e-5;

% We save last three system solutions
% is the error estimation
y_ = zeros(length(initVals),3);
dt_ = zeros(1,3);
y_(:,end) = initVals;
dt_(end) = dt_optimal_;

sysIndex_ = 0;

if multirate
    ysize = erkSize;
    yTypical = P_ERK.solver.yTypical;
else
    ysize = length(yTypical);
end

while t_ < tEnd
    
    sysIndex_ = sysIndex_ + 1;
    
    if(rem(sysIndex_, 1000)==0) % for displaying progress
        
        sysIndex_, toc, t_
        
    end
    
    % Calculate the solution
    [sol_,  PERS] = solve_sys([t_ t_ + dt_optimal_],...
        relTol, step_rejected_, PERS);
    
    [dt_optimal_, step_rejected_, P_ERK.controller] = ec_h211b(...
        [y_(1:ysize,:) sol_(1:ysize)], dt_, sysIndex_,...
        relTol, yTypical, PERS.CELL.controller.eEst, P_ERK.controller);

    if (step_rejected_)
        dt_(end) = dt_optimal_;
        sysIndex_ = sysIndex_ - 1;
    else
        t_ = t_ + dt_(end);
        dt_=circshift(dt_,[0,-1]);
        y_=circshift(y_,[0,-1]);
        dt_(end) = dt_optimal_;
        y_(:,end) = sol_;
    end
    
    % Ensure that we dont leave the interval
    if t_ + dt_optimal_ > tEnd  
        dt_optimal_ = tEnd-t_;
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
