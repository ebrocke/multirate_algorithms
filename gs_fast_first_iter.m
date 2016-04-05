function SYSTEM = gs_fast_first_iter(t, relTol, SYSTEM)

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

t_ = t(1);
ii_ = 0;
% H_ is a macro time step
cellfirst=1;
erkfirst=0;
while t_ < t(2)
    
    [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
        @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK);
    [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
        @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL);
    
    
    h_cell = get_h_optimal(SYSTEM.CELL);
    h_erk = get_h_optimal(SYSTEM.ERK);
    
    if(rem(ii_, 1000)==0) % for displaying progress
        [toc, t_, h_cell, h_erk ]
    end
    
    % the component with the smallest h_ is solved first over H_
    if (h_cell <= h_erk)
        
        % full rollback if component order has been changed
        if (erkfirst && step_rejected)
            SYSTEM.ERK = update([t_ t_], SYSTEM.ERK);
        end
        
        H_ = h_erk;
        t_ = SYSTEM.ERK.controller.t(end);
        
        SYSTEM.CELL.sys.isolver_hdl = @solve_adaptive_step;
        SYSTEM.ERK.sys.isolver_hdl = @solve_one_step;
        
        [out SYSTEM] = cell_first([t_ t_+H_],...
            {t_erk, y_erk},...
            relTol, SYSTEM);
        
        cellfirst = 1;
        erkfirst = 0;
    else
        % full rollback if component order has been changed
        if (cellfirst && step_rejected)
            SYSTEM.CELL = update([t_ t_], SYSTEM.CELL);
        end
        
        H_ = h_cell;
        t_ = SYSTEM.CELL.controller.t(end);
        
        SYSTEM.ERK.sys.isolver_hdl = @solve_adaptive_step;
        SYSTEM.CELL.sys.isolver_hdl = @solve_one_step;
        
        [out SYSTEM] = erk_first([t_ t_+H_],...
            {t_cell, y_cell}, ...
            relTol, SYSTEM);
        
        cellfirst = 0;
        erkfirst = 1;
    end
    
    if (any(isnan(out)))
        step_rejected = true;
    else
        step_rejected = false;
    end
    
    ii_ = ii_ +1;
    
end

end