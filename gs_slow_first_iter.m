function SYSTEM = gs_slow_first_iter(t, relTol, SYSTEM)
step_rejected = false;
t_ = t(1);
ii_ = 0;
% H_ is a macro time step

while t_ < t(2)
    
    % retrieve last 3 values of the exchanged variables
    [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
        @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK);
    [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
        @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL);
    
    h_cell = get_h_optimal(SYSTEM.CELL);
    h_erk = get_h_optimal(SYSTEM.ERK);
    
    if(rem(ii_, 1000)==0) % for displaying progress
        [toc, t_, h_cell h_erk ]
    end
    
    if (h_cell > h_erk)
        H_ = h_cell;
        t_ = SYSTEM.CELL.controller.t(end);
        SYSTEM.CELL.sys.isolver_hdl = @solve_one_step;
        SYSTEM.ERK.sys.isolver_hdl = @solve_adaptive_step;
        
        [out SYSTEM] = CellFirst([t_ t_+H_],...
            {t_erk, y_erk, t_cell, y_cell},...
            relTol, step_rejected, SYSTEM);
    else
        H_ = h_erk;
        t_ = SYSTEM.ERK.controller.t(end);
        SYSTEM.CELL.sys.isolver_hdl = @solve_adaptive_step;
        SYSTEM.ERK.sys.isolver_hdl = @solve_one_step;
        
        [out SYSTEM] = ERKFirst([t_ t_+H_],...
            {t_erk, y_erk, t_cell, y_cell},...
            relTol, step_rejected, SYSTEM);
    end
    
    if (any(isnan(out)))
        step_rejected = true;
    else
        step_rejected = false;
    end
    
    ii_ = ii_ +1;
end
end