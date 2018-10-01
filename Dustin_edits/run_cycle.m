function [Weight,param] = run_cycle(options,mission)
A1 = std_atmosphere(options.height);%Ambient conditions as a function of altitude
[OTM,A2,A3,A4,A5,O1,O2,O3,O4,O5] = OxygenModule(options,A1);
[FC,E1] = oxy_fuelcell(options,O5);
[HL,F1,F2,F3,F4,E2,E3,E4] = HeatLoop(options,FC,OTM,E1);
[Weight,param] = NetParam(options,FC,OTM,HL);

param.states = {'A1',A1;'A2',A2;'A3',A3;'A4',A4;'A5',A5;'E1',E1;'E2',E2;'E3',E3;'E4',E4;'F1',F1;'F2',F2;'F3',F3;'F4',F4;'O1',O1;'O2',O2;'O3',O3;'O4',O4;'O5',O5;};
end