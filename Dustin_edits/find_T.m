function T = find_T(Flow, H)
Flow_2 = Flow;
Flow_3 = Flow;
error = 1;
count = 0;
while max(max(abs(error)))>1e-3
    count = count+1;
    Flow_2.T = Flow.T+.01;
%     Flow_3.T = Flow.T-1;
    H_flow = enthalpy(Flow); %property(Flow,'h','kJ');
    dT_dH = .01./(enthalpy(Flow_2) - H_flow); %(property(Flow_2,'h','kJ') - H_flow);
%     dT_dH = 0.5*(.01./(property(Flow_2,'h','kJ') - H_flow) + 1./(H_flow - property(Flow_3,'h','kJ')));
    error = dT_dH.*(H -  H_flow);
    error(isinf(error)) = 0;
    Flow.T = Flow.T + error;
    if count>25
        Flow.T = Flow_3.T;
        error = 0;
    end
end
T = Flow.T;
end