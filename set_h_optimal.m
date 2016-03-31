function SYSTEM =set_h_optimal(h, SYSTEM)
ii_ = SYSTEM.stats.acceptedIter + 1;
SYSTEM.sol.dt(ii_) = h;
end