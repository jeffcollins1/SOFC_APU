function param = run_cycle(options,mission)
A1.T = -0.0065*options.height + 14.987 + 273.1; %Ambient Temperature Kelvin;
A1.P = 107*exp(-0.0001*options.height)- 10; %Ambient Pressure as a function of altitude
A1.O2 = .21*options.airflow;
A1.N2 = .79*options.airflow;
[OTM,A1,A2,A3,A4,A5,O1,O2,O3] = OxygenModule(options,A1);
[FC,E1] = FuelCell(options,O3);
[F1,F2,F3,F4,F5,E2,E3,E4,E5] = HeatLoop(options,FC,OTM,E1);
param = NetParam(options,mission,FC,OTM);


param.states = {'A1',A1;'A2',A2;'A3',A3;};

