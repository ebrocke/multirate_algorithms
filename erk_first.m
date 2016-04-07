function [out SYSTEM] = erk_first(t, vars, relTol, SYSTEM, ec_fn)
ec_erk_local = ec_fn{1};
ec_cell_local = ec_fn{2};

    function [H STEP_REJECTED PERSISTENT] = ec_cell( DT, E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_cell_local(...
            DT, E_EST, PERSISTENT, SYSTEM.ERK);
    end

    function [H STEP_REJECTED PERSISTENT] = ec_erk( DT,  E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_erk_local(...
            DT, E_EST, PERSISTENT, SYSTEM.CELL);
    end
 

SYSTEM.CELL.controller.fn = @ec_cell;
SYSTEM.ERK.controller.fn = @ec_erk;

t_erk = vars{1};
y_erk = vars{2};
t_cell = vars{3};
y_cell = vars{4};
isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

[Y1, SYSTEM.ERK] = isolver_erk(t,...
    {t_cell, y_cell}, ...
    relTol, ...
    SYSTEM.ERK);

if any(isnan(Y1))
    erk_vars = [NaN NaN];
else
    [fracTilde, caiTilde] = feval(SYSTEM.ERK.sys.exch_hdl, Y1);
    erk_vars = [fracTilde caiTilde*1e3];
end

[Y2, SYSTEM.CELL] = isolver_cell(t,...
    {[-(t(2)-t(1)), t_erk], [erk_vars; y_erk]}, ...
    relTol,...
    SYSTEM.CELL);
out = [Y1; Y2];
end