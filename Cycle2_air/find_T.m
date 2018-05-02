function T = find_T(Flow, H)
Flow_2 = Flow;
Flow_2.T = Flow.T+1;
dT_dH = 1/(property(Flow_2,'h','kJ')-property(Flow,'h','kJ'));

error = 1;
while abs(error)>1e-3
    error = (H - property(Flow,'h','kJ'))/H;
%     if abs(dT_dH*error*H)>50
%         Flow_2 = Flow;
%         Flow_2.T = Flow.T+1;
%         dT_dH = .5*dT_dH + .5/(property(Flow_2,'h','kJ')-property(Flow,'h','kJ'));
%     end
    Flow.T = Flow.T + 0.75*dT_dH*error*H;
end
T = Flow.T;
end