function sim_combine2j2()
clear global; close all;
%dbstop if warning; dbstop if error;
%%%% Generating map structure
addpath(genpath('methods')); %,genpath('kkit'), genpath('common'));
global NGATES
global MODE CONST 
global iterMethod

MODE=3;
CONST = load('modelconst.mat');

%%%% Initiate solvers
organization = ...              % Organization for solve decoupled systems
    {'Jac','GSwErkfirst','GSwCellfirst'};
iterMethod = organization{3};

%%%% Chosing solver parameters
multirate = true;
% Erk system
solvErkNM   = @BDF2_DEF;%@BDF2_DEF;%@BDF2_AS;%@RK4; 
isolvErk = @solve_one_step;
% Cell system
solvCellNM   = @BDF2_DEF;%@BDF2_DEF;%@BDF2_DEF;%@BDF2_AS;%@CrankNicStagg; 
isolvCell = @solve_fixed_step;

odeslv = @solver_adaptive;
%odeslv = @ode_solver_fixed;%@ode_solver_h211b; %ode_solver_fixed
slv_param = 1e-4;   % Relative tolerance for adaptive [ @ode_solver_h211b (1e-3)]
                    % or number of steps for fixed [ @ode_solver_fixed (1e5)]
%slv_param = 2e5; 
                    
%%%% Chosing model parameters
T_SIM = 2;         % Total run time [s]

%%%% Set up the parametrs
setup_cell_parameters();
setup_erk_parameters();

%%%% Initial sate variables
cell_init_vals = init_cell_variables(true);
%ca_pump=false when simulated model_v.3
ca_pump = true;
erk_init_vals = init_erk_variables(ca_pump, 0);
erk_size = length(erk_init_vals);

%%%% Interface between two systems
[~, erk_exch]= solve_erk_handle();
[~, cell_exch] = solve_cell_handle();

%%%% Run the simulation
tic;
[t1, y_erk, t2, y_cell, stats_erk, stats_cell] = odeslv([0 T_SIM], ...
    [erk_init_vals cell_init_vals], slv_param, ...
    erk_size, ...
    {solvErkNM erk_exch isolvErk}, ...
    {solvCellNM cell_exch isolvCell},...
    multirate);
runTime = toc

%%%% Saving the data

mkdir(sprintf('../Data/%s',date));
%%%% Define strings for saved filename
if multirate
    m_ = '_mr';
else
    m_ = '_sr';
end
filename = sprintf('../Data/%s/%s.%s.%g.%s%s.%g.%g',date,...
    func2str(solvErkNM), func2str(solvCellNM), MODE, iterMethod, m_, slv_param, T_SIM);
fn1 = sprintf('%s.mat',filename);
fn_y1 = sprintf('%s_erk_ca_pmapk_ka.mat',filename);
fn_y2 = sprintf('%s_hh_gates.mat', filename);
fn_y3 = sprintf('%s_hh_v.mat',filename);
%save the workspace except y vector
save(fn1,'-regexp','^(?!(y_erk|y_cell)$).')
% save biochemical values
erk_vars = y_erk([1,6,9],:);
save(fn_y1,'erk_vars');
% save gate probabilities
hh_gates = y_cell(1:NGATES,:);
save(fn_y2,'hh_gates');
% save voltages
hh_v = y_cell((1+NGATES):end,:);
save(fn_y3,'hh_v');
%
%save(filename);
end

% function [frac, cai] = get_exch_dummy(t) % in SI
% global EVARS
%     frac = interp1(EVARS.TTILDE,EVARS.FRACTILDE,t);
%     cai = interp1(EVARS.TTILDE,EVARS.CAITILDE,t);
% end

function [caflux] = get_exch_dummy(t) % in SI
global EVARS
    caflux = interp1(EVARS.TTILDE,EVARS.CAFLUXTILDE,t);
end