function [tint, y]=ode_solver_fixed(solve_sys,tint, ...
    init_vals, Nsys, erkSize, erkSys, cellSys)
%global MODE
%global Nsys Nhh Nerk
%global y_vec y_vec_erk
global dt sysIndex stats%hh_iter erk_iter ode_cell_iter ode_erk_iter

% y_vec = zeros(length(init_vals(erk_size+1:end)),Nsys*Nhh);
% y_vec(:,1) = init_vals(erk_size+1:end);
% y_vec_erk = zeros(length(init_vals(1:erk_size)),Nsys*Nerk);
% y_vec_erk(:,1) = init_vals(1:erk_size);
% 
% hh_iter = zeros(1,Nsys*Nhh);
% erk_iter = zeros(1,Nsys*Nerk);
% ode_cell_iter = 0;
% ode_erk_iter  = 0;
step_rejected = false;


yTypical = importdata('yTypicalSolution.txt',';'); % Typical solution

% Persistent now handled in a more concise way

P_ERK.init = true;

P_ERK.yTypical = yTypical(1:erkSize);

P_CELL.init = true;

P_CELL.yTypical = yTypical(erkSize+1:end);

PERS.ERK = P_ERK;

PERS.CELL = P_CELL;


%Defining time intervals and step sizes
t0 = tint(1); T_SIM = tint(2);

tint = linspace(t0, T_SIM, Nsys);

% we have a fixed step size
% the history of 3 values is enough
dt(1:3) = tint(2)-tint(1);
% Set up initial values etc.
y = zeros(length(init_vals),1);
y(:,1) = init_vals;


for sysIndex = 1:Nsys-1
    
    if(rem(sysIndex,1000)==0) % for displaying progress
        sysIndex
        toc
        sysIndex*dt(end)
    end
    
   %     [frac_tilde, cai_tilde]  = feval(erk_sys{2},y(1:erk_size,sys_index));
   %     erk_tilde               = [frac_tilde cai_tilde];
    % For the first steps or MODE=1, use constant predictors
%     if (sys_index < 2) || (sys_index < 3 && MODE==3) || (MODE==1)  
%         [frac_tilde, cai_tilde]  = feval(erk_sys{2},y(1:erk_size,sys_index));
%         erk_tilde               = [frac_tilde cai_tilde*1e3];
%         
%         [ca_flux_tilde] = feval(cell_sys{2},y(erk_size+1:end,sys_index),cai_tilde);
%         cell_tilde      = ca_flux_tilde;
%     elseif (MODE == 2) %linear polynomial predictor
%         erk_tilde = zeros(2,2);
%         cell_tilde = zeros(2,1);
%         for jj = [1 0]
%             [frac_tilde, cai_tilde] = feval(erk_sys{2},y(1:erk_size,sys_index-jj));
%             erk_tilde(end-jj,:) = [frac_tilde cai_tilde*1e3];
%             
%             [ca_flux_tilde]    = feval(cell_sys{2},y(erk_size+1:end,sys_index-jj),cai_tilde);
%             cell_tilde(end-jj) = ca_flux_tilde;
%         end    
%     elseif (MODE == 3)% second order polynomial predictor
%         erk_tilde = zeros(3,2);
%         cell_tilde = zeros(3,1);
%         for jj = [2 1 0]
%             [frac_tilde, cai_tilde] = feval(erk_sys{2},y(1:erk_size,sys_index-jj));
%             erk_tilde(end-jj,:) = [frac_tilde cai_tilde*1e3];
%             
%             [ca_flux_tilde]    = feval(cell_sys{2},y(erk_size+1:end,sys_index-jj),cai_tilde);
%             cell_tilde(end-jj) = ca_flux_tilde;
%         end
%     else
%         disp('You need to specify the MODE {1 ,2, 3}');
%     end
%     states = y(:,sys_index);
    [out, stat, PERS]=solve_sys([tint(sysIndex) tint(sysIndex+1)], y, erkSize,...
        erkSys, cellSys, 1e-5, PERS, step_rejected);
    if any(isnan(out))
        disp(['Slowly convergent Newton method, h = ', num2str(dt(end))])
        break;
    end
    stats = stats+stat;
    y(:,sysIndex+1) = out;
end
% if(~sum(erk_iter))
%     erk_iter = Nsys*Nerk;
% end
% if(~sum(hh_iter))
%     hh_iter = Nsys*Nhh;
% end

fprintf('%-20s \t %-8s \t %-8s\n',sprintf('System'), 'Erk','HH')
fprintf('%-20s \t %-8d \t %-8d\n','numjac calls', stats(1,1), stats(2,1))
fprintf('%-20s \t %-8d \t %-8d\n','ode calls by numjac', stats(1,2), stats(2,2))
fprintf('%-20s \t %-8d \t %-8d\n','ode calls by iterations', stats(1,3), stats(2,3))

% fprintf('%-20s \t %-8s \t %-8s\n',sprintf('System'), 'Erk','HH')
% fprintf('%-20s \t %-8d \t %-8d\n','Solver calls',sum(erk_iter),sum(hh_iter))
% fprintf('%-20s \t %-8d \t %-8d\n','Failed solver calls',...
%     sum(erk_iter)-Nerk*Nsys,sum(hh_iter)-Nsys*Nhh)
% fprintf('%-20s \t %-8d \t %-8d\n','ode calls',sum(ode_erk_iter),sum(ode_cell_iter))
% 
end