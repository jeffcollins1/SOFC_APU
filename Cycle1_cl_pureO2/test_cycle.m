%%  Test Cycle
%Create options to test
options.OTM_area = 1e4; %membrane area in m^2
options.SOFC_area = 1e4; %SOFC active area





Param = run_cycle(options,mission);
