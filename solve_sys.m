function [out, SYSTEM] = solve_sys(t, relTol, step_rejected, SYSTEM)

global iterMethod

% H = t(2)-t(1);
% isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
% isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

    function v = get_erk_exchaged_vector(erk_sol, ~)
        [frac cai] = feval(SYSTEM.ERK.sys.exch_hdl,erk_sol);
        v = [frac cai*1e3]; %mM
    end



% t0 - time back from current (t=0) at which the variables are requested
    function v = get_cell_exchaged_vector(cell_sol, t0)
        CAI = approximate(t_erk,y_erk(:,2),t0);
        [caFlux] = feval(SYSTEM.CELL.sys.exch_hdl, cell_sol, CAI*1e-3);
        v = [0 caFlux];
    end


% single rate integration
if(strcmp(iterMethod,'Jac'))
    isolver_cell = @solve_one_step;
    isolver_erk = @solve_one_step;
    t_ = t(1);
    ii_ = 0;
    % H_ is a macro time step
    while t_ < t(2)
        [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
            @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK);
        [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
            @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL);
        
        if(rem(ii_, 1000)==0) % for displaying progress
            [toc, t_, get_h_optimal(SYSTEM.ERK) ]
        end    
        
        % ERK component calculates H_=h_
        H_ = get_h_optimal(SYSTEM.ERK);
        t_ = SYSTEM.ERK.controller.t(end);
        
        
        SYSTEM.CELL = set_h_optimal(H_,SYSTEM.CELL);
        SYSTEM.CELL.controller.eEst  = 0;
        
        erkTilde = approximate(t_erk, y_erk, -H_);
        [Y2, SYSTEM.CELL] = isolver_cell([t_ t_+H_], ...
            {[-H_ t_erk], [erkTilde; y_erk]}, ...
            relTol,...
            SYSTEM.CELL);
       
        
        SYSTEM.ERK.controller.eEst = SYSTEM.CELL.controller.eEst;
        cellTilde = approximate(t_cell, y_cell, -H_ );
        [Y1, SYSTEM.ERK] = isolver_erk([t_ t_+H_],...
            {[-H_ t_cell], [cellTilde; y_cell]},...
            relTol,...
            SYSTEM.ERK);
        
        if (any(isnan(Y1)) && ~any(isnan(Y2)))
            SYSTEM.CELL = update([t_ t_],SYSTEM.CELL);
        end
        
        out = [Y1; Y2];
        if (any(isnan(out)))
            step_rejected = true;
        else
            step_rejected = false;
        end
        
        ii_ = ii_ +1;
    end

% slowest component integrated first
% (predicted time step is the largest)
elseif(strcmp(iterMethod,'GSSlowFirst'))
    t_ = t(1);
    while t_ < t(2)
        ii_ = max(SYSTEM.ERK.stats.acceptedIter, ...
            SYSTEM.CELL.stats.acceptedIter);
        if(rem(ii_, 1000)==0) % for displaying progress
            [toc, t_, SYSTEM.CELL.controller.h SYSTEM.ERK.controller.h ]
        end
        
        % retrieve last 3 values of the exchanged variables
        [y_erk t_erk SYSTEM.ERK] = get_exchanged_vector(...
            @get_erk_exchaged_vector, step_rejected, SYSTEM.ERK);
        [y_cell t_cell SYSTEM.CELL] = get_exchanged_vector( ...
            @get_cell_exchaged_vector, step_rejected, SYSTEM.CELL);
        
        if (SYSTEM.CELL.controller.h > SYSTEM.ERK.controller.h)
            H_ = SYSTEM.ERK.controller.h;
            [out SYSTEM] = CellFirst([t_ t_+H_],...
                {t_erk, y_erk, t_cell, y_cell},...
                relTol, step_rejected, SYSTEM);
        else
            H_ = SYSTEM.CELL.controller.h;
            [out SYSTEM] = ERKFirst([t_ t_+H_],...
                {t_erk, y_erk, t_cell, y_cell},...
                relTol, step_rejected, SYSTEM);
        end
        t_ = t_ + H_;
    end
    % fastest component integrated first
    % (predicted time step is the smallest)
elseif(strcmp(iterMethod,'GSFastFirst'))
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

            [out SYSTEM] = CellFirst([t_ t_+H_],...
                {t_erk, y_erk, t_cell, y_cell},...
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
            [out SYSTEM] = ERKFirst([t_ t_+H_],...
                {t_erk, y_erk, t_cell, y_cell}, ...
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
    
else
    fprintf('Wrong iteration method');
end
end

function [out SYSTEM] = ERKFirst(t, vars, relTol, SYSTEM)
t_erk = vars{1};
y_erk = vars{2};
t_cell = vars{3};
y_cell = vars{4};
isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

%cellTilde = approximate(t_cell, y_cell, -H);

%     [Y1, SYSTEM.ERK] = isolver_erk(t,...
%         {[-H t_cell], [cellTilde; y_cell]}, ...
%         relTol, step_rejected,...
%         SYSTEM.ERK);
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
SYSTEM.CELL.controller.eEst = SYSTEM.ERK.controller.eEst;
[Y2, SYSTEM.CELL] = isolver_cell(t,...
    {[-(t(2)-t(1)) t_erk], [erk_vars; y_erk]}, ...
    relTol,...
    SYSTEM.CELL);
out = [Y1; Y2];
end

function [out SYSTEM] = CellFirst(t, vars, relTol, SYSTEM)
t_erk = vars{1};
y_erk = vars{2};
%t_cell = vars{3};
%y_cell = vars{4};
isolver_cell = SYSTEM.CELL.sys.isolver_hdl;
isolver_erk = SYSTEM.ERK.sys.isolver_hdl;

% if we use extrapolation as first
% and then interpolation for each micro time step,
% then we need to recalculate the whole step
% in case the step is rejected
% Thus we use extrapolation for each micro time step

%     [Y2, SYSTEM.CELL] = isolver_cell(t, ...
%         {[-H t_erk], [erkTilde; y_erk]}, ....
%         relTol, step_rejected, ...
%         SYSTEM.CELL);
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
SYSTEM.ERK.controller.eEst = SYSTEM.CELL.controller.eEst;
[Y1, SYSTEM.ERK] = isolver_erk( t, ...
    { [], cell_vars},...
    relTol, ...
    SYSTEM.ERK);

out = [Y1; Y2];
end

