function v = get_erk_exchaged_vector(SYSTEM, erk_sol, ~, ~)
[frac cai] = feval(SYSTEM.sys.exch_hdl,erk_sol);
v = [frac cai*1e3]; %mM
end

