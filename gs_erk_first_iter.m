function SYSTEM = gs_erk_first_iter(t, relTol, SYSTEM)

    function [H STEP_REJECTED PERSISTENT] = ec_cell( DT, E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_h211b(...
            DT, max(E_EST, SYSTEM.ERK.controller.eEst), PERSISTENT);
    end

    function [H STEP_REJECTED PERSISTENT] = ec_erk( ~, ~, PERSISTENT)
        STEP_REJECTED = false;
        H = 0; % setup manually for erk component later on
    end


step_rejected = false;
isolver_cell = @solve_one_step;
isolver_erk = @solve_one_step;
SYSTEM.CELL.controller.fn = @ec_cell;
SYSTEM.ERK.controller.fn = @ec_erk;

t_ = t(1);
ii_ = 0;

while t_ < t(2)
    
    [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
        @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK, {});
    [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
        @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL, {t_erk, y_erk});
    
    % CELL component defines macro time step H_=h_
    H_ = get_h_optimal(SYSTEM.CELL);
    % CELL component defines current time
    t_ = SYSTEM.CELL.controller.t(end);
    % CELL controller calculates optimal time step,
    % thus we have to set h manually
    SYSTEM.ERK = set_h_optimal(H_, SYSTEM.ERK);
    
    if(rem(ii_, 1000)==0) % for displaying progress
        [toc, t_, H_ ]
    end
    
    cellTilde = approximate(t_cell, y_cell, -H_ );
        
%     [out, SYSTEM] = erk_first([t_ t_+H_], ...
%         {[-H_ t_cell], [cellTilde; y_cell]},...
%         relTol, SYSTEM);
    
    [Y1, SYSTEM.ERK] = isolver_erk( [t_ t_+H_], ...
        {[-H_ t_cell], [cellTilde; y_cell]},...
        relTol, ...
        SYSTEM.ERK);
    
    if any(isnan(Y1))
        erk_vars = [NaN NaN];
        %display('CellFirst:: varsTilde is NaN')
    else
        [fracTilde, caiTilde] = feval(SYSTEM.ERK.sys.exch_hdl, Y1);
        erk_vars = [fracTilde caiTilde*1e3];
    end
    
    [Y2, SYSTEM.CELL] = isolver_cell([t_ t_+H_],...
        {[], erk_vars}, ...
        relTol,...
        SYSTEM.CELL);
    
    out = [Y1; Y2];
    
    if (any(isnan(out)))
        step_rejected = true;
        SYSTEM.ERK = update([t_ t_],SYSTEM.ERK);
    else
        step_rejected = false;
    end
    
    ii_ = ii_ +1;
    
end
end