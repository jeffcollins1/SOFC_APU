function T = find_T(Flow, H)
Flow_2 = Flow;
Flow_2.T = Flow.T+1;
[n1,n2] = size(Flow.T)
dT_dH = 1./((enthalpy(Flow_2) - enthalpy(Flow)));%(property(Flow_2,'h','kJ')-property(Flow,'h','kJ'));

error = 1;
while max(max(abs(error)))>1e-3
    error = (H - enthalpy(Flow))./H; %property(Flow,'h','kJ'))./H;
    Flow.T = Flow.T + 0.75*dT_dH.*error.*H;
end
T = Flow.T;
end