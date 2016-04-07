function [out SYSTEM] = cell_first(t, vars, relTol, SYSTEM)

    function [H STEP_REJECTED PERSISTENT] = ec_cell( DT, E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_h211b(...
            DT, max(E_EST, SYSTEM.ERK.controller.eEst), PERSISTENT);
    end

    function [H STEP_REJECTED PERSISTENT] = ec_erk( DT,  E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_h211b(...
            DT, max(E_EST, SYSTEM.CELL.controller.eEst), PERSISTENT);
    end

SYSTEM.CELL.controller.fn = @ec_cell;
SYSTEM.ERK.controller.fn = @ec_erk;

t_erk = vars{1};
y_erk = vars{2};
t_cell = vars{3};
y_cell = vars{4};
isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

[Y2, SYSTEM.CELL] = isolver_cell(t, ...
    {t_erk, y_erk}, ....
    relTol, ...
    SYSTEM.CELL);

if any(isnan(Y2))
    cell_vars = [NaN NaN];
else
    erkTilde = approximate(t_erk, y_erk, -(t(2)-t(1)));
    [caFlux] = feval(SYSTEM.CELL.sys.exch_hdl,...
        Y2, erkTilde(2)*1e-3);
    cell_vars = [0 caFlux];
end
    
[Y1, SYSTEM.ERK] = isolver_erk( t, ...
    { [-(t(2)-t(1)), t_cell], [cell_vars; y_cell]},...
    relTol, ...
    SYSTEM.ERK);

out = [Y1; Y2];
end