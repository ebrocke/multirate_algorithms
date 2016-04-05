function [t_erk, y_erk, t_cell, y_cell, stats_erk, stats_cell]=...
    solver_adaptive(t, ...
    initVals, ...
    relTol, ...
    erkSize, ...
    erkSys, cellSys,...
    multirate)

global MODE iterMethod


%tEnd = t(2);     % Ending time
%t_ = t(1);         % Current time

%step_rejected = false;

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
P_ERK.solver.step_rejected = false;
P_ERK.solver.yTypical = yTypical(1:erkSize);
ss = length(P_ERK.solver.yTypical);
P_ERK.solver.fac_ = zeros(ss,1);   
P_ERK.solver.delta_old_ = zeros(ss,3);
P_ERK.solver.new_jac_ = true;
P_ERK.solver.j_= zeros(ss,ss);

% controller working variables
P_ERK.controller.eEst = 0;
P_ERK.controller.init = 0;
%P_ERK.controller.h = 1e-5;
P_ERK.controller.t = t(1);
P_ERK.controller.eEst_ = 0;
P_ERK.controller.eEstVec_ =  [NaN NaN NaN]; 
P_ERK.controller.rhofac_ = 0;
P_ERK.controller.eEstVec = [NaN 0.9 0.9];
P_ERK.controller.rhofac = 0;

% system working variables
P_ERK.sys.t_exch = zeros(1,MODE); %store last exchanged values
P_ERK.sys.y_exch = zeros(MODE,2); % in a reversed order of time
% handles that provides flexibility in solving the system
P_ERK.sys.method_hdl = erkSys{1}; % BDF_DEF BDF_AS RK4 CrankNicStagg; 
P_ERK.sys.exch_hdl = erkSys{2}; % system function that provides exchanged variables
%P_ERK.sys.isolver_hdl = erkSys{3};% how to solve interval (used in multirate)
P_ERK.sys.ode_hdl = @ode_erk; % handle to the system ode functions
% store system solution
P_ERK.sol.y = zeros(erkSize,mem_size);
P_ERK.sol.dt = zeros(1,mem_size);
P_ERK.sol.dt(1) = 1e-5;
P_ERK.sol.y(:,1) = initVals(1:erkSize)';


P_CELL.solver.init = true;
P_CELL.solver.yTypical = yTypical(erkSize+1:end);
P_CELL.solver.step_rejected = false;
ss = length(P_CELL.solver.yTypical);
P_CELL.solver.fac_ = zeros(ss,1);   
P_CELL.solver.delta_old_ = zeros(ss,3);
P_CELL.solver.new_jac_ = true;
P_CELL.solver.j_= zeros(ss,ss);

P_CELL.stats.acceptedIter = 0;
P_CELL.stats.refinedIter = 0;
P_CELL.stats.rejectedIter = 0;
P_CELL.stats.numjac = 0;
P_CELL.stats.n_ode_numjac = 0;
P_CELL.stats.n_ode_iter = 0;
P_CELL.sol.y = zeros(length(initVals)-erkSize,mem_size);
P_CELL.sol.dt = zeros(1,mem_size);
P_CELL.sol.dt(1) = 1e-5;
P_CELL.sol.y(:,1) = initVals(erkSize+1:end)';
P_CELL.controller.init = 0;
P_CELL.controller.eEst = 0;
%P_CELL.controller.h = 1e-5;
P_CELL.controller.t = t(1);
P_CELL.controller.eEst_ = 0;
P_CELL.controller.eEstVec_ = [NaN NaN NaN]; 
P_CELL.controller.rhofac_ = 0;
P_CELL.controller.eEstVec = [NaN 0.9 0.9];
P_CELL.controller.rhofac = 0;

P_CELL.sys.method_hdl = cellSys{1};
P_CELL.sys.exch_hdl = cellSys{2};
%P_CELL.sys.isolver_hdl = cellSys{3};
P_CELL.sys.ode_hdl = @ode_cell;
P_CELL.sys.t_exch = zeros(1,MODE);
P_CELL.sys.y_exch = zeros(MODE,2);
P_CELL.sys.h = 1e-5;
P_CELL.sys.m = 1;
P_CELL.sys.eEst  = 0;

SYSTEM.ERK = P_ERK;
SYSTEM.CELL = P_CELL;
SYSTEM.yTypical = yTypical;

%yTypical = abs(yTypical);
%H_ = 1e-5;

% We save last three system solutions
% is the error estimation
%y_ = zeros(length(initVals),3);
%dt_ = zeros(1,3);
%y_(:,end) = initVals;
%dt_(end) = H_;


% if multirate
%    ysize = erkSize;
%    yTypical = P_ERK.solver.yTypical;
% else
%     ysize = length(yTypical);
% end


if(strcmp(iterMethod,'Jac'))
    SYSTEM = jac_iter(t, relTol, SYSTEM); %single rate
elseif(strcmp(iterMethod,'GSCELLFirst'))
    SYSTEM = gs_cell_first_iter(t, relTol, SYSTEM); %single rate
elseif(strcmp(iterMethod,'GSERKFirst'))
    SYSTEM = gs_erk_first_iter(t, relTol, SYSTEM); %single rate
elseif(strcmp(iterMethod,'GSSlowFirst'))
    SYSTEM = gs_slow_first_iter(t, relTol, SYSTEM); %multirate
elseif(strcmp(iterMethod,'GSFastFirst'))
    SYSTEM = gs_fast_first_iter(t, relTol, SYSTEM); %multirate
else
    fprintf('Wrong iteration method');
end


% fill in return values
i_ = SYSTEM.ERK.stats.acceptedIter;
y_erk = SYSTEM.ERK.sol.y(:, 1:i_+1);
t_erk = [0 cumsum(SYSTEM.ERK.sol.dt(1:i_))];
stats_erk = SYSTEM.ERK.stats;

i_ = SYSTEM.CELL.stats.acceptedIter;
y_cell = SYSTEM.CELL.sol.y(:, 1:i_+1);
t_cell = [0 cumsum(SYSTEM.CELL.sol.dt(1:i_))];
stats_cell = SYSTEM.CELL.stats;

% 
 fprintf('%-20s \t %-8s \t %-8s\n',sprintf('System'), 'Erk','HH')
% 
 fprintf('%-20s \t %-8d \t %-8d\n','numjac calls', ...
     SYSTEM.ERK.stats.numjac, SYSTEM.CELL.stats.numjac)
% 
 fprintf('%-20s \t %-8d \t %-8d\n','n ode by numjac', ...
     SYSTEM.ERK.stats.n_ode_numjac, SYSTEM.CELL.stats.n_ode_numjac)
% 
 fprintf('%-20s \t %-8d \t %-8d\n','n ode by iterations',...
     SYSTEM.ERK.stats.n_ode_iter, SYSTEM.CELL.stats.n_ode_iter)
%
  fprintf('%-20s \t %-8d \t %-8d\n','accepted calls', ...
     SYSTEM.ERK.stats.acceptedIter, SYSTEM.CELL.stats.acceptedIter)
%
  fprintf('%-20s \t %-8d \t %-8d\n','refined calls', ...
     SYSTEM.ERK.stats.refinedIter, SYSTEM.CELL.stats.refinedIter)
 %
  fprintf('%-20s \t %-8d \t %-8d\n','rejected calls', ...
     SYSTEM.ERK.stats.rejectedIter, SYSTEM.CELL.stats.rejectedIter)
end
