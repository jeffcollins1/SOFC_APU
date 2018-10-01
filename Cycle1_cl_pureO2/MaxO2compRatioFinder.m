O3.T = 323*ones(10,10); %Compressor inlet temperature
O3.P = linspace(30,50,10); 
O3.O2 = 0.21*ones(10,1)*ones(1,10);
O3.N2 = 0.79*ones(10,10); 
options.P_fc = 1000*ones(10,10);
options.C2_eff = 0.80*ones(10,10); 



[O4,C2_work] = compressor(O3,options.P_fc,options.C2_eff)

Pratio = O4.P./O3.P; 

Final = [O4.T(1,1:10);Pratio(1,1:10)]