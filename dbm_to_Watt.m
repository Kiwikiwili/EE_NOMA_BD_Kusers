function [P_watt] = dbm_to_Watt(P_dbm)

P_watt=1e-3*10.^(P_dbm/10);

end

