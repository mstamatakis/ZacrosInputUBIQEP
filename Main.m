clear all
clc

global Rgas SurfSpecsQT EquivGasSurfSpecsNames...
    EqGasSurfaceSpecsFormH gammaRh BondIndex A_s_Forw A_s_Back T T0...
    beta_Forw beta_Back MW

Rgas = 8.314/4.184/1000; %[kcal/mol/K]
T0 = 300; %[K]

T = 500 + 273.15; %[K]
deltaT = T - T0; %[K]

gammaRh = 2.49e-8; %kmol/m2

MW = [1.00794,17.00734,18.01528,28.0101,44.0095,12.0107,13.01864,14.02658,15.03452,45.01744,16.04246,15.9994,2.01588,18.01528,28.0101,44.0095]/1000; %[kg/mol]
EquivGasSurfSpecsNames = {'H*','OH*','H2O*','CO*','CO2*','C*','CH*','CH2*','CH3*','COOH*','CH4','O*','H2','H2O','CO','CO2'};


for i = 1:length(EquivGasSurfSpecsNames)
    EqGasSurfaceSpecsFormH(i) = EquivGasFormationEnthalpy(EquivGasSurfSpecsNames{i},T);
end

SurfSpecsQT0_ZeroCoverage = [62.3 70 10.8 38.5 5.2 159 151.2 109.3 42.4 62.2 6 100]; %[kcal/mol]

IE = zeros (12,5);
IE(1,1) = -2.5;
IE(1,4) = -3.7;
IE(2,2) = -25;
IE(2,3) = 25;
IE(2,5) = -33;
IE(3,2) = 25;
IE(3,3) = -4.5;
IE(4,1) = -3.7;
IE(4,4) = -15;
IE(12,5) = -26;

teta_H = 0.0;
teta_OH = 0.0;
teta_H2O = 0;
teta_CO = 0;
teta_O = 0;

teta = [teta_H , teta_OH , teta_H2O , teta_CO , teta_O];

f = [1.5 , 2 , 2.5 , 2 , 2 , 1.5 , 2 , 2.5 , 2.5 , 2.5 , 2 , 1.5];

SurfSpecsQT = SurfSpecsQT0_ZeroCoverage + (IE * teta')' - f * Rgas * deltaT;

irxn = 0;

irxn = irxn + 1;
GasReact{irxn} = {'H2'};
SurfReact{irxn} = {};
GasProd{irxn} = {};
SurfProd{irxn} = {'H*','H*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 7.73e-1;
A_s_Back(irxn) = 5.56e11;
beta_Forw(irxn) = 0.9387;
beta_Back(irxn) = -0.4347;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'H*','OH*'};
BondIndex(irxn) = 0.55;
A_s_Forw(irxn) = 1.15e11;
A_s_Back(irxn) = 3.6e8;
beta_Forw(irxn) = 0.0281;
beta_Back(irxn) = 1.2972;

irxn = irxn + 1;
GasReact{irxn} = {'H2O'};
SurfReact{irxn} = {};
GasProd{irxn} = {};
SurfProd{irxn} = {'H2O*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 7.72e-2;
A_s_Back(irxn) = 2.06e13;
beta_Forw(irxn) = 1.4067;
beta_Back(irxn) = -1.8613;

irxn = irxn + 1;
GasReact{irxn} = {'CO'};
SurfReact{irxn} = {};
GasProd{irxn} = {};
SurfProd{irxn} = {'CO*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 5.00e-1;
A_s_Back(irxn) = 5.65e12;
beta_Forw(irxn) = -2.000;
beta_Back(irxn) = 1.9879;

irxn = irxn + 1;
GasReact{irxn} = {'CO2'};
SurfReact{irxn} = {};
GasProd{irxn} = {};
SurfProd{irxn} = {'CO2*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 3.67e-1;
A_s_Back(irxn) = 7.54e10;
beta_Forw(irxn) = -2.3294;
beta_Back(irxn) = 2.1831;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CO2*','H*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CO*','OH*'};
BondIndex(irxn) = 0.7;
A_s_Forw(irxn) = 6.42e-1;
A_s_Back(irxn) = 5.64e-1;
beta_Forw(irxn) = 0.030;
beta_Back(irxn) = 0.03006;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'COOH*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CO*','OH*'};
BondIndex(irxn)= 0.5;
A_s_Forw(irxn) = 1.07e12;
A_s_Back(irxn) = 9.37e11;
beta_Forw(irxn) = -0.4123;
beta_Back(irxn) = 0.4123;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'COOH*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CO2*','H*'};
BondIndex(irxn)= 2.2;
A_s_Forw(irxn) = 1.0e10;
A_s_Back(irxn) = 9.99e9;
beta_Forw(irxn) = -0.4424;
beta_Back(irxn) = 0.4424;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CO*','H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'COOH*','H*'};
BondIndex(irxn) = 0.9;
A_s_Forw(irxn) = 3.34e9;
A_s_Back(irxn) = 1.2e7;
beta_Forw(irxn) = -0.2222;
beta_Back(irxn) = 0.2223;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CO2*','H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'COOH*','OH*'};
BondIndex(irxn) = 0.01;
A_s_Forw(irxn) = 1.90e13;
A_s_Back(irxn) = 5.60e10;
beta_Forw(irxn) = -0.1922;
beta_Back(irxn) = 0.1922;

irxn = irxn + 1;
GasReact{irxn} = {'CH4'};
SurfReact{irxn} = {};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH3*','H*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 5.72e-1;
A_s_Back(irxn) = 7.72e10;
beta_Forw(irxn) = 0.7883;
beta_Back(irxn) = -0.7883;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH3*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH2*','H*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 2.49e10;
A_s_Back(irxn) = 2.57e9;
beta_Forw(irxn) = 0.0862;
beta_Back(irxn) = -0.0862;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH2*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH*','H*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 5.5e10;
A_s_Back(irxn) = 7.27e9;
beta_Forw(irxn) = -0.1312;
beta_Back(irxn) = 0.1312;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'C*','H*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 4.58e12;
A_s_Back(irxn) = 2.18e11;
beta_Forw(irxn) = -0.2464;
beta_Back(irxn) = 0.2464;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH3*','O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH2*','OH*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 2.96e11;
A_s_Back(irxn) = 3.38e10;
beta_Forw(irxn) = -0.1906;
beta_Back(irxn) = 0.1906;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH2*','H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH3*','OH*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 5.73e10;
A_s_Back(irxn) = 1.74e9;
beta_Forw(irxn) = -0.7208;
beta_Back(irxn) = 0.7208;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CH*','H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH2*','OH*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 6.49e11;
A_s_Back(irxn) = 1.541e10;
beta_Forw(irxn) = -0.5033;
beta_Back(irxn) = 0.5033;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'C*','H2O*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'CH*','OH*'};
BondIndex(irxn) = 0.5;
A_s_Forw(irxn) = 9.74e11;
A_s_Back(irxn) = 6.41e10;
beta_Forw(irxn) = -0.3882;
beta_Back(irxn) = 0.3882;

irxn = irxn + 1;
GasReact{irxn} = {};
SurfReact{irxn} = {'CO*','H*'};
GasProd{irxn} = {};
SurfProd{irxn} = {'C*','OH*'};
BondIndex(irxn) = 0.45;
A_s_Forw(irxn) = 1.18e12;
A_s_Back(irxn) = 7.60e12;
beta_Forw(irxn) = 0.2944;
beta_Back(irxn) = -0.2944;

CRLF = [char(13) char(10)];
fileID = fopen('fileout.txt','w');
fileZacrosID = fopen('mechanism_input.dat','w');
fprintf (fileZacrosID,['mechanism' CRLF]);
fprintf (fileZacrosID,[CRLF]);

for i=1:irxn
    ReactionNumber(i) = i;
    UBIType(i) = IdReaction(GasReact{i}, SurfReact{i}, GasProd{i}, SurfProd{i});
    y{i} = EactUBIQEP(GasReact{i}, SurfReact{i}, GasProd{i}, SurfProd{i}, UBIType(i),ReactionNumber(i));
    EactForw{i} = y{i}(1);
    EactBack{i} = y{i}(2);
    kForw{i}= y{i}(3);
    kBack{i}= y{i}(4);
    
    disp (['Reaction ' , num2str(ReactionNumber(i))]);
    disp ([GasReact{i}, SurfReact{i}, '-->', GasProd{i}, SurfProd{i}]);
    disp (['Reaction_Type: ' , num2str(UBIType(i))]);
    disp (['Eact_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(EactForw{i}) ,' kcal/mol']);
    disp (['Eact_React', num2str(ReactionNumber(i)), '_Back: ', num2str(EactBack{i}) ,' kcal/mol']);
    disp (['k_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(kForw{i},'%1.2e')]);
    disp (['k_React', num2str(ReactionNumber(i)), '_Back: ', num2str(kBack{i},'%1.2e')]);
    disp ('**************************************************************');
    
    fprintf (fileID,['Reaction ' , num2str(ReactionNumber(i)) CRLF],'\n');
    
    if length(GasReact{i})==1 && length(SurfProd{i})==1
        fprintf (fileID,[char(cell2mat(GasReact{i}(1))),' + *']);
    elseif length(GasReact{i})==1 && length(SurfProd{i})==2
        fprintf (fileID,[char(cell2mat(GasReact{i}(1))),' + 2*']);
    elseif length(GasReact{i})==2
        fprintf (fileID,[char(cell2mat(GasReact{i}(1))),' + ',char(cell2mat(GasReact{i}(2)))]);
    end
    
    if length(SurfReact{i})==1
        fprintf (fileID,[char(cell2mat(SurfReact{i}(1))),' + *']);
    elseif length(SurfReact{i})==2
        fprintf (fileID,[char(cell2mat(SurfReact{i}(1))),' + ',char(cell2mat(SurfReact{i}(2)))]);
    end
    
    fprintf (fileID,[' --> ']);
    
    if length(GasProd{i})==1
        fprintf (fileID,[char(cell2mat(GasProd{i}(1))) CRLF]);
    elseif length(GasProd{i})==2
        fprintf (fileID,[char(cell2mat(GasProd{i}(1))),' + ',char(cell2mat(GasProd{i}(2))) CRLF]);
    end
    
    if length(SurfProd{i})==1
        fprintf (fileID,[char(cell2mat(SurfProd{i}(1))) CRLF]);
    elseif length(SurfProd{i})==2
        fprintf (fileID,[char(cell2mat(SurfProd{i}(1))),' + ',char(cell2mat(SurfProd{i}(2))), CRLF]);
    end
    
    fprintf (fileID,['Reaction_Type: ' , num2str(UBIType(i)) CRLF]);
    fprintf (fileID,['Eact_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(EactForw{i}) ,' kcal/mol' CRLF]);
    fprintf (fileID,['Eact_React', num2str(ReactionNumber(i)), '_Back: ', num2str(EactBack{i}) ,' kcal/mol' CRLF]);
    fprintf (fileID,['k_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(kForw{i},'%1.2e') CRLF]);
    fprintf (fileID,['k_React', num2str(ReactionNumber(i)), '_Back: ', num2str(kBack{i},'%1.2e') CRLF]);
    fprintf (fileID,['**************************************************************' CRLF]);
    fprintf (fileID,[' ' CRLF]);

    fprintf (fileZacrosID,['reversible_step Reaction' , num2str(i) CRLF]);
    fprintf (fileZacrosID,['   sites ' , num2str(max(length(SurfReact{i}))) CRLF]);

    fprintf (fileZacrosID,['end_reversible_step' CRLF]);
    fprintf (fileZacrosID,['' CRLF]);
    
    
    
    
    
    
end

fprintf (fileZacrosID,['end_mechanism' CRLF]);

fclose(fileID);
fclose(fileZacrosID);

