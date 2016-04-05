% t0 - time back from current (t=0) at which the variables are requested
function v = get_cell_exchaged_vector(SYSTEM, cell_sol, t0, erk_sol)
t_erk = erk_sol{1};
y_erk = erk_sol{2};
CAI = approximate(t_erk,y_erk(:,2),t0);
[caFlux] = feval(SYSTEM.sys.exch_hdl, cell_sol, CAI*1e-3);
v = [0 caFlux];
end