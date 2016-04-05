function [out SYSTEM] = erk_first(t, vars, relTol, SYSTEM)
%t_erk = vars{1};
%y_erk = vars{2};
t_cell = vars{1};
y_cell = vars{2};
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
    {[], erk_vars}, ...
    relTol,...
    SYSTEM.CELL);
out = [Y1; Y2];
end