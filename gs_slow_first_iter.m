function SYSTEM = gs_slow_first_iter(t, relTol, SYSTEM)
 eEst = [NaN NaN NaN];
 
    function [H STEP_REJECTED PERSISTENT] = ec( DT, E_EST, PERSISTENT, R_SYSTEM)
        eEst = circshift(eEst,[0 1]);
        eEst(1) = R_SYSTEM.controller.eEst;
        if(length(DT) > 2)
            eEst_ = approximate([0 DT(end-1)],eEst(1:2)',-DT(end));
        else 
            eEst_ = eEst(1);
        end
        [H,  STEP_REJECTED, PERSISTENT] = ec_h211b(...
            DT,  max(E_EST, eEst_), PERSISTENT);
        if STEP_REJECTED
            eEst = circshift(eEst,[0 -1]);
        end
    end

step_rejected = false;
t_ = t(1);
ii_ = 0;

%SYSTEM.CELL.controller.fn = @ec_cell;
%SYSTEM.ERK.controller.fn = @ec_erk;
% H_ is a macro time step

while t_ < t(2)
    
    % retrieve last 3 values of the exchanged variables
    [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
        @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK, {});
    [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
        @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL, {t_erk, y_erk});
    
    h_cell = get_h_optimal(SYSTEM.CELL);
    h_erk = get_h_optimal(SYSTEM.ERK);
    
    if(rem(ii_, 100)==0) % for displaying progress
        [toc, t_, h_cell, h_erk ]
    end
    
    if (h_cell > h_erk)
        H_ = h_cell;
        t_ = SYSTEM.CELL.controller.t(end);
        SYSTEM.CELL.sys.isolver_hdl = @solve_one_step;
        SYSTEM.ERK.sys.isolver_hdl = @solve_adaptive_step;
        %erkTilde = approximate(t_erk, y_erk, -H_);
        [out SYSTEM] = cell_first([t_ t_+H_],...
            {t_erk, y_erk, t_cell, y_cell},...
            relTol, SYSTEM, {@ec, @ec});
    else
        H_ = h_erk;
        t_ = SYSTEM.ERK.controller.t(end);
        SYSTEM.CELL.sys.isolver_hdl = @solve_adaptive_step;
        SYSTEM.ERK.sys.isolver_hdl = @solve_one_step;
        %cellTilde = approximate(t_cell, y_cell, -H_ );
        [out SYSTEM] = erk_first([t_ t_+H_],...
            {t_erk, y_erk, t_cell, y_cell},...
            relTol, SYSTEM, {@ec, @ec});
    end
    
    if (any(isnan(out)))
        step_rejected = true;
    else
        step_rejected = false;
    end
    
    ii_ = ii_ +1;
end
end