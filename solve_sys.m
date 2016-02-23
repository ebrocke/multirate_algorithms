function [out, SYSTEM] = solve_sys(t, relTol, step_rejected, SYSTEM)

global iterMethod

H = t(2)-t(1);
isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

    function v = get_erk_exchaged_vector(erk_sol, ~)
        [frac cai] = feval(SYSTEM.ERK.sys.exch_hdl,erk_sol);
        v = [frac cai*1e3]; %mM
    end

[y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
    @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK);

% t0 - time back from current (t=0) at which the variables are requested
    function v = get_cell_exchaged_vector(cell_sol, t0)
        CAI = approximate(t_erk,y_erk(:,2),t0);
        [caFlux] = feval(SYSTEM.CELL.sys.exch_hdl, cell_sol, CAI*1e-3);
        v = [0 caFlux];
    end

[y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
    @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL);

if(strcmp(iterMethod,'Jac'))
    cellTilde = approximate(t_cell, y_cell, -H );
    [Y1, SYSTEM.ERK] = isolver_erk(t,...
        {[-H t_cell], [cellTilde; y_cell]},...
        relTol, step_rejected,...
        SYSTEM.ERK);
    
    
    erkTilde = approximate(t_erk, y_erk, -H);
    [Y2, SYSTEM.CELL] = isolver_cell(t, ...
        {[-H t_erk], [erkTilde; y_erk]}, ...
        relTol, step_rejected,...
        SYSTEM.CELL);
    
    out = [Y1; Y2];
    
elseif(strcmp(iterMethod,'GSwErkfirst'))
    cellTilde = approximate(t_cell, y_cell, -H);
    
    [Y1, SYSTEM.ERK] = isolver_erk(t,...
        {[-H t_cell], [cellTilde; y_cell]}, ...
        relTol, step_rejected,...
        SYSTEM.ERK);
    
    if any(isnan(Y1))
        erk_vars = [NaN NaN];
    else
        [fracTilde, caiTilde] = feval(SYSTEM.ERK.sys.exch_hdl, Y1);
        erk_vars = [fracTilde caiTilde*1e3];
    end
    
    [Y2, SYSTEM.CELL] = isolver_cell(t,...
        {[-H t_erk], [erk_vars; y_erk]}, ...
        relTol, step_rejected,...
        SYSTEM.CELL);
    
    out = [Y1; Y2];
    
elseif(strcmp(iterMethod,'GSwCellfirst'))
    
    % if we use extrapolation as first
    % and then interpolation for each micro time step,
    % then we need to recalculate the whole step 
    % in case the step is rejected
    % Thus we use extrapolation for each micro time step
    
%     [Y2, SYSTEM.CELL] = isolver_cell(t, ...
%         {[-H t_erk], [erkTilde; y_erk]}, ....
%         relTol, step_rejected, ...
%         SYSTEM.CELL);
    [Y2, SYSTEM.CELL] = isolver_cell(t, ...
        {t_erk, y_erk}, ....
        relTol, step_rejected, ...
        SYSTEM.CELL);
    
    if any(isnan(Y2))
        cell_vars = [NaN NaN];
    else
        erkTilde = approximate(t_erk, y_erk, -H);
        [caFlux] = feval(SYSTEM.CELL.sys.exch_hdl,...
            Y2, erkTilde(2)*1e-3);
        cell_vars = [0 caFlux];
    end
    
    [Y1, SYSTEM.ERK] = isolver_erk( t, ...
        {[], cell_vars},...
        relTol, step_rejected, ...
        SYSTEM.ERK);
    
    out = [Y1; Y2];
else
    fprintf('Wrong iteration method');
end
end