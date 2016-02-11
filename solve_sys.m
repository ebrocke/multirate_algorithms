function [out, SYSTEM] = solve_sys(t, relTol, SYSTEM, step_rejected)

global iterMethod 
H = t(2)-t(1);

interval_solver_cell = @solve_interval_one_step; %@solve_interval_adaptive_step; %@solve_interval_fixed_step
interval_solver_erk = @solve_interval_one_step; 

    function v = get_erk_exchaged_vector(erk_sol, ~)
        [frac cai] = feval(SYSTEM.ERK.sys.exch_hdl,erk_sol);
        v = [frac cai*1e3]; %mM
    end

fun_hdl = @get_erk_exchaged_vector;
[y_erk dt_erk SYSTEM.ERK] = get_exchanged_vector(SYSTEM.ERK, step_rejected,fun_hdl);

    % t0 - time back from current (t=0) at which the variables are requested
    function v = get_cell_exchaged_vector(cell_sol, t0)
%         if (cum_dt > 0 && abs(cum_dt - dt_erk(end)) > 1e-15) 
%             % if multirate we need to interpolate CAI values
%             % for intermendiete dt
%             if cum_dt < dt_erk(end)
%                 CAI = interpolate(y_erk(:,2),dt_erk,dt_erk(end)-cum_dt);
%             else
%                 CAI = interpolate(flipud(y_erk(:,2)), fliplr(dt_erk), cum_dt-dt_erk(end)) ;
%             end
%         else
%             % we need to pick up CAI value corresponding to the dt
%             cum = [0 cumsum(fliplr(dt_erk))]-cum_dt;
%             I = find(abs(cum)<1e-10);
%             if (length(I) > 1 || length(I) == 0)
%                 point = 0;
%             end
%            
%             CAI = y_erk(end-I+1,2);
%         end
        CAI = approximate(dt_erk,y_erk(:,2),t0);
        [caFlux] = feval(SYSTEM.CELL.sys.exch_hdl, cell_sol, CAI*1e-3);
        v = [ 0 caFlux];
    end

fun_hdl = @get_cell_exchaged_vector;
[y_cell dt_cell SYSTEM.CELL] = get_exchanged_vector(SYSTEM.CELL, step_rejected,fun_hdl);

if(strcmp(iterMethod,'Jac'))
    cellTilde = extrapolate(y_cell, [dt_cell H] );
    [Y1, SYSTEM.ERK] = interval_solver_erk(SYSTEM.ERK.sys.slv_hdl, ...
        t, @ode_erk, {cellTilde, [dt_cell H]}, relTol, SYSTEM.ERK, step_rejected);
    
   
    erkTilde = extrapolate(y_erk, [dt_erk H]);
    [Y2, SYSTEM.CELL] = interval_solver_cell(SYSTEM.CELL.sys.slv_hdl,...
        t, @ode_cell, {erkTilde, [dt_erk H]}, relTol, SYSTEM.CELL, step_rejected);
    
    out = [Y1; Y2];
    
elseif(strcmp(iterMethod,'GSwErkfirst'))
    cellTilde = extrapolate(y_cell, [dt_cell H]);
    [Y1,SYSTEM.ERK] = interval_solver_erk(SYSTEM.ERK.sys.slv_hdl, ...
        t, @ode_erk, {cellTilde, [dt_cell H]}, relTol, SYSTEM.ERK, step_rejected);
    
    if any(isnan(Y1))
        erkTilde = [NaN NaN];
    else
        %the last values in erkTilde are known
        [fracTilde, caiTilde] = feval(SYSTEM.ERK.sys.exch_hdl, Y1);
        %erkTilde should have the history of previous values
        % (these values are used for interpolation in the multirate
        % integration)
        erkTilde = circshift(y_erk,-1);
        erkTilde(end,:)   = [fracTilde caiTilde*1e3];
        
    end
    [Y2, SYSTEM.CELL] = interval_solver_cell(SYSTEM.CELL.sys.slv_hdl,...
        t, @ode_cell, {erkTilde, [dt_erk H]}, relTol, SYSTEM.CELL, step_rejected);
    
    out = [Y1; Y2];
    
elseif(strcmp(iterMethod,'GSwCellfirst'))
  
    erkTilde = approximate(dt_erk, y_erk, -H);
    
    [Y2, SYSTEM.CELL] = interval_solver_cell(SYSTEM.CELL.sys.slv_hdl,...
        t, @ode_cell, {erkTilde, [dt_erk H]}, relTol , SYSTEM.CELL, step_rejected);
    
    if any(isnan(Y2))
        caFluxTilde = NaN;
    else
        [caFluxTilde] = feval(SYSTEM.CELL.sys.exch_hdl, Y2, erkTilde(end,2)*1e-3);
    end
    
    [Y1, SYSTEM.ERK] = interval_solver_erk(SYSTEM.ERK.sys.slv_hdl, ...
            t, @ode_erk, {[0 caFluxTilde], []}, relTol, SYSTEM.ERK, step_rejected);
    out = [Y1; Y2];
else
    fprintf('Wrong iteration method');
end
end