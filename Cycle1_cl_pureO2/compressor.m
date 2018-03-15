function [A2,C1_work] = compressor(A1,options.P_non_perm,options.C1_eff)
[P,T] = intake(options.height,options.velocity)
O2.P = P;
O2.T = T;
N2.P = P;
N2.T = T;
hin = 0.21*enthalpy(O2) + 0.79*enthalpy(N2);
sin = 0.21*entropy(O2) + 0.79*entropy(N2);
A2.P = options.P_non_perm;

