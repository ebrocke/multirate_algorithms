function h =get_h_optimal(SYSTEM)
ii_ = SYSTEM.stats.acceptedIter + 1;
h = SYSTEM.sol.dt(ii_);
end