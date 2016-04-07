function SYSTEM = jac_iter(t, relTol, SYSTEM)

    function [H STEP_REJECTED PERSISTENT] = ec_cell( ~, ~, PERSISTENT)
        STEP_REJECTED = false;
        H = 0; % setup manually for cell component later on
    end

    function [H STEP_REJECTED PERSISTENT] = ec_erk( DT, E_EST, PERSISTENT)
        [H,  STEP_REJECTED, PERSISTENT] = ec_h211b(...
    DT, max(E_EST, SYSTEM.CELL.controller.eEst), PERSISTENT);
    end

step_rejected = false;
isolver_cell = @solve_one_step;
isolver_erk = @solve_one_step;
SYSTEM.CELL.controller.fn = @ec_cell;
SYSTEM.ERK.controller.fn = @ec_erk;

t_ = t(1);
ii_ = 0;

% H_ is a macro time step
while t_ < t(2)
    [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
        @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK, {});
    [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
        @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL, {t_erk,y_erk});
    
    % ERK component defines macro time step H_=h_
    H_ = get_h_optimal(SYSTEM.ERK);
    % ERK component defines current time
    t_ = SYSTEM.ERK.controller.t(end);
    % ERK controller calculates optimal time step,
    % thus we have to set h manually
    SYSTEM.CELL = set_h_optimal(H_,SYSTEM.CELL);
    
    if(rem(ii_, 1000)==0) % for displaying progress
        [toc, t_, H_ ]
    end
   
    %erkTilde = approximate(t_erk, y_erk, -H_);
    %{[-H_ t_erk], [erkTilde; y_erk]},
    [Y2, SYSTEM.CELL] = isolver_cell([t_ t_+H_], ...
        {t_erk, y_erk}, ...
        relTol,...
        SYSTEM.CELL);
    

    %cellTilde = approximate(t_cell, y_cell, -H_ );
    %{[-H_ t_cell], [cellTilde; y_cell]},...
    [Y1, SYSTEM.ERK] = isolver_erk([t_ t_+H_],...
        {t_cell, y_cell},...
        relTol,...
        SYSTEM.ERK);
    
    out = [Y1; Y2];
    
    if (any(isnan(out)))
        step_rejected = true;
        SYSTEM.CELL = update([t_ t_],SYSTEM.CELL);
    else
        step_rejected = false;
    end
    
    ii_ = ii_ +1;
end
end
