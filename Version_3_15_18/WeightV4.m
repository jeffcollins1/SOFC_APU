
%Weight Analysis
function [Pden,WeightHX,WeightCompressor,WeightTurbine,WeightFC,WeightOTM,WeightSystem] = WeightV4(QE1_2,QO1_2,QE4_F3,QE2_3,WC2,Preq,PT1,WC1,Cells,Membrane);
if QO1_2 >0
Qtransfer = QE1_2+ QO1_2 + QE4_F3 + QE2_3; % kW
PdenHX = 15; % kW/kg
WeightHX = Qtransfer/PdenHX; % kg
Pdenmotor = 20; %kW/kg
WeightMotor = Preq/Pdenmotor; 
Pdenturbine = 8; %kW/kg
WeightTurbine = PT1/Pdenturbine; %kg
Wdencompressor = 5; %kW/kg
WeightMembrane = .0473; %Weight of one 10cm x 10cm OTM membrane and housing
WcompTotal = WC1 + WC2; 
WeightOTM = Membrane*WeightMembrane;
WeightCompressor = WcompTotal/Wdencompressor; %kg
WeightCell = 0.06534; %Weight of Cell, 2 Thermiculite seals and 1 flow plate,
SWvessel = 0.3; %kg per liter
WeightFC = Cells*(WeightCell);
WeightH2insulator = 270; %kg
WeightSystem = WeightMotor + WeightCompressor + WeightTurbine + WeightFC +  WeightHX + WeightOTM;
MassFuel = 1660; %kg
% MaxPowerRunTime = MassFuel/(mFlowH2new*3600); %Length of time available at max power, hours
Pden = Preq/WeightSystem;
else 
    Pden = 0;
    WeightHX = 0;
    WeightCompressor =0;
    WeightTurbine = 0;
    WeightFC = 0;
    WeightSystem = 1;
end
% WeightVec(z) = SystemWeight(z);
% EffSystemVec(z) = EffSystem(z);
% a = max(PdenSystemVec)-0.01;
% find(PdenSystemVec>a,z);
% WeightHX = WeightHX/WeightSystem(z); 
% WeightFC = WeightFC/WeightSystem(z);
% Weight = SystemWeight(z);
end