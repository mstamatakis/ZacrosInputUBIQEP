clear all
close all
fclose all
clc

global Rgas SurfSpecsQT EquivGasSurfSpecsNames...
    EqGasSurfaceSpecsFormH gammaRh BondIndex A_s_Forw A_s_Back T T0...
    beta_Forw beta_Back MW

random_seed = 71543 ;

Rgas = 8.314/4.184/1000; %[kcal/mol/K]
T0 = 300; %[K]

T = 500 + 273.15; %[K]
deltaT = T - T0; %[K]

gammaRh = 2.49e-8; %kmol/m2

MW = [1.00794,17.00734,18.01528,28.0101,44.0095,12.0107,13.01864,14.02658,15.03452,45.01744,16.04246,15.9994,2.01588,18.01528,28.0101,44.0095]/1000; %[kg/mol]
EquivGasSurfSpecsNames = {'H*','OH*','H2O*','CO*','CO2*','C*','CH*','CH2*','CH3*','COOH*','CH4','O*','H2','H2O','CO','CO2'};
GasSpecies = {'CH4','H2','H2O','CO','CO2'};
SurfaceSpecies = {'H*','OH*','H2O*','CO*','CO2*','C*','CH*','CH2*','CH3*','COOH*','O*'};

P = 1 ; %bar %random value
Molec_Fract_CH4 = 0.2; %random value
Molec_Fract_H2 = 0.2; %random value
Molec_Fract_H2O = 0.3; %random value
Molec_Fract_CO = 0.2; %random value
Molec_Fract_CO2 = 0.1; %random value
Molecular_Fractions = [Molec_Fract_CH4 Molec_Fract_H2 Molec_Fract_H2O Molec_Fract_CO Molec_Fract_CO2];
Partial_Pressure = P*Molecular_Fractions; %bar

for i = 1:length(EquivGasSurfSpecsNames)
    EqGasSurfaceSpecsFormH(i) = EquivGasFormationEnthalpy(EquivGasSurfSpecsNames{i},T);
end

Species_for_Q = {'H*','OH*','H2O*','CO*','CO2*','C*','CH*','CH2*','CH3*','COOH*', 'CH4' ,'O*'};
Species_for_Q_for_energetics = {'H','OH','H2O','CO','CO2','C','CH','CH2','CH3','COOH', 'CH4' ,'O'};

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
NameReaction{irxn} = {'H2-dissociative-ads',};
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
NameReaction{irxn} = {'H2O*-dissociation'};
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
NameReaction{irxn} = {'H2O-ads',};
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
NameReaction{irxn} = {'CO-ads'};
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
NameReaction{irxn} = {'CO2-ads'};
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
NameReaction{irxn} = {'CO2*-deoxidation'};
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
NameReaction{irxn} = {'COOH*-dissociation-v1'};
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
NameReaction{irxn} = {'COOH*-dissociation-v2'};
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
NameReaction{irxn} = {'CO*-hydroxylation'};
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
NameReaction{irxn} = {'CO2*-hydroxylation'};
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
NameReaction{irxn} = {'CH4-dissociative-ads'};
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
NameReaction{irxn} = {'CH3*-dehydrogenation-v1'};
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
NameReaction{irxn} = {'CH2*-dehydrogenation'};
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
NameReaction{irxn} = {'CH*-dehydrogenation'};
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
NameReaction{irxn} = {'CH3*-dehydrogenation-v2'};
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
NameReaction{irxn} = {'CH2*-hydrogenation'};
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
NameReaction{irxn} = {'CH*-hydrogenation'};
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
NameReaction{irxn} = {'C*-hydrogenation',};
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
NameReaction{irxn} = {'CO-deoxidation'};
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

fileID = fopen('summary_reactions.txt','w');

fileZacrosID_Mechanism = fopen('mechanism_input.dat','w');
fprintf (fileZacrosID_Mechanism,['mechanism' CRLF]);
fprintf (fileZacrosID_Mechanism,[CRLF]);
fprintf (fileZacrosID_Mechanism,['#########################################################' CRLF]);

for i=1:irxn
    ReactionNumber(i) = i;
    UBIType(i) = IdReaction(GasReact{i}, SurfReact{i}, GasProd{i}, SurfProd{i});
    y{i} = EactUBIQEP(GasReact{i}, SurfReact{i}, GasProd{i}, SurfProd{i}, UBIType(i),ReactionNumber(i));
    EactForw{i} = y{i}(1);
    EactBack{i} = y{i}(2);
    PreExponForw{i} = y{i}(3);
    PreExponBack{i} = y{i}(4);
    kForw{i}= y{i}(5);
    kBack{i}= y{i}(6);
    
    disp (['Reaction ' , num2str(ReactionNumber(i))]);
    disp ([GasReact{i}, SurfReact{i}, '-->', GasProd{i}, SurfProd{i}]);
    disp (['Reaction_Type: ' , num2str(UBIType(i))]);
    disp (['Eact_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(EactForw{i}) ,' kcal/mol']);
    disp (['Eact_React', num2str(ReactionNumber(i)), '_Back: ', num2str(EactBack{i}) ,' kcal/mol']);
    disp (['Pre_Expon', num2str(ReactionNumber(i)), '_Forw: ', num2str(PreExponForw{i},'%1.2e')]);
    disp (['Pre_Expon', num2str(ReactionNumber(i)), '_Back: ', num2str(PreExponBack{i},'%1.2e')]);
    disp (['k_React', num2str(ReactionNumber(i)), '_Forw: ', num2str(kForw{i},'%1.2e')]);
    disp (['k_React', num2str(ReactionNumber(i)), '_Back: ', num2str(kBack{i},'%1.2e')]);
    disp ('#########################################################');
    
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
    fprintf (fileID,['#########################################################' CRLF]);
    fprintf (fileID,[' ' CRLF]);
    
    
    %From here we create the mechanism input for Zacros
    fprintf (fileZacrosID_Mechanism,['' CRLF]);
    fprintf (fileZacrosID_Mechanism,['reversible_step ', cell2mat(NameReaction{i}) CRLF]);
    
    if ~isempty (GasReact{i})
        fprintf (fileZacrosID_Mechanism,['  gas_reacs_prods  ', cell2mat(GasReact{i}),' -' , num2str(length(GasReact{i})) CRLF]);
    end
    fprintf (fileZacrosID_Mechanism,['  sites ', num2str(max([length(GasReact{i}),length(SurfReact{i}),length(GasProd{i}),length(SurfProd{i})])) CRLF]);
    if max([length(GasReact{i}),length(SurfReact{i}),length(GasProd{i}),length(SurfProd{i})])==2
        fprintf (fileZacrosID_Mechanism,['  neighboring 1-2' CRLF]);
    end
    
    fprintf (fileZacrosID_Mechanism,['  initial' CRLF]);
    if max([length(GasReact{i}),length(SurfReact{i}),length(GasProd{i}),length(SurfProd{i})])==1
        if length(GasReact{i})==1
            fprintf (fileZacrosID_Mechanism,['    1 * 1' CRLF]); %adsorption
        end
    else
        if isempty(SurfReact{i})
            for ind=1:length(SurfProd{i})
                fprintf (fileZacrosID_Mechanism,['    ', num2str(ind), ' * 1'  CRLF]); %dissociative adsorption
            end
        elseif length(SurfReact{i})==1
            
            fprintf (fileZacrosID_Mechanism,['    1 ',cell2mat(SurfReact{i}), ' 1'  CRLF]);
            fprintf (fileZacrosID_Mechanism,['    2 *    1'  CRLF]); %Surface Reaction of 1 surface species
            
            
        elseif length(SurfReact{i})==2
            for ind=1:length(SurfReact{i})
                fprintf (fileZacrosID_Mechanism,['    ', num2str(ind), ' ',cell2mat(SurfReact{i}(ind)), ' 1'  CRLF]); %Surface Reaction of 2 surface species
            end
            
            
            
        end
    end
    
    fprintf (fileZacrosID_Mechanism,['  final' CRLF]);
    for ind=1:length(GasProd{i})
        fprintf (fileZacrosID_Mechanism,['    ', num2str(ind), ' ', cell2mat(GasProd{i}(ind)), ' 1' CRLF]);
    end
    for ind=1:length(SurfProd{i})
        fprintf (fileZacrosID_Mechanism,['    ', num2str(ind), ' ', cell2mat(SurfProd{i}(ind)), ' 1'  CRLF]);
    end
    
    fprintf (fileZacrosID_Mechanism,['  site_types ']);
    for indexsite = 1:max([length(GasReact{i}),length(SurfReact{i}),length(GasProd{i}),length(SurfProd{i})])
        fprintf (fileZacrosID_Mechanism,['StTp1 ']);
    end
    fprintf (fileZacrosID_Mechanism,[CRLF]);
    
    fprintf (fileZacrosID_Mechanism,['  pre_expon ', num2str(PreExponForw{i},'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  pe_ratio ', num2str(PreExponForw{i}/PreExponBack{i},'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  activ_eng ', num2str(EactForw{i},'%1.2e')  CRLF]);
    
    if ~isempty(GasReact{i})
        if length(SurfProd{i})==1
            fprintf (fileZacrosID_Mechanism,['  prox_factor ', sprintf('%1.1f',0.0) CRLF]);
        elseif length(SurfProd{i})==2
            fprintf (fileZacrosID_Mechanism,['  prox_factor ', sprintf('%1.1f',1.0) CRLF]);
        end
    else
        fprintf (fileZacrosID_Mechanism,['  prox_factor ', sprintf('%1.1f',0.5) CRLF]);
    end
    
    fprintf (fileZacrosID_Mechanism,['end_reversible_step' CRLF]);
    fprintf (fileZacrosID_Mechanism,['' CRLF]);
    fprintf (fileZacrosID_Mechanism,['#########################################################' CRLF]);
end

fprintf (fileZacrosID_Mechanism,['' CRLF]);
fprintf (fileZacrosID_Mechanism,['end_mechanism' CRLF]);

fclose(fileID);
fclose(fileZacrosID_Mechanism);

%From here we create the simulation_input file
fileZacrosID_Simulation = fopen('simulation_input.dat','w');

fprintf (fileZacrosID_Simulation,['random_seed               ', num2str(random_seed) , ' # random value' CRLF]);
fprintf (fileZacrosID_Simulation,['' CRLF]);
fprintf (fileZacrosID_Simulation,['temperature               ', sprintf('%1.2f',T) CRLF]);
fprintf (fileZacrosID_Simulation,['pressure                  ', sprintf('%1.1f',P), ' # random value' CRLF]);
fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['n_gas_species             ', num2str(length(GasSpecies)) CRLF]);

fprintf (fileZacrosID_Simulation,['gas_specs_names             ']);
for indexgasspecies = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,[cell2mat(GasSpecies(indexgasspecies)),'     ']);
end
fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['gas_energies              ']);

%I choose as reference species CH4 H2 H2O

Reference_Gas_Species = {'CH4' , 'H2' , 'H2O'};

for indexRefGasSpecies = 1:length(Reference_Gas_Species)
    Matrix_H_ref_gas_species(indexRefGasSpecies) = EqGasSurfaceSpecsFormH(SpecIndx(Reference_Gas_Species{indexRefGasSpecies},EquivGasSurfSpecsNames));
end

MatrixCompRef = [1 0 0; 4 2 2; 0 0 1]; % in rows: C H O , in column: CH4, H2, H2O
vett_comp_CH4 = [0 0 0]'; %if ref species [0 0 0] if not components respect to [C H O]
vett_comp_H2 = [0 0 0]'; %if ref species [0 0 0] if not components respect to [C H O]
vett_comp_H2O = [0 0 0]'; %if ref species [0 0 0] if not components respect to [C H O]
vett_comp_CO = [1 0 1]'; %if ref species [0 0 0] if not components respect to [C H O]
vett_comp_CO2 = [1 0 2]'; %if ref species [0 0 0] if not components respect to [C H O]
vett_comp = [vett_comp_CH4 vett_comp_H2 vett_comp_H2O vett_comp_CO vett_comp_CO2];
for indexcomp = 1:length(vett_comp)
    coeff_form_energy(:,indexcomp) = inv(MatrixCompRef) * vett_comp(:,indexcomp);
end

for indexformenergy = 1:length(GasSpecies)
    if indexformenergy <= length(Reference_Gas_Species)
        FormEnergy(indexformenergy) = 0;
    else
        FormEnergy(indexformenergy) = (EqGasSurfaceSpecsFormH(SpecIndx(GasSpecies{indexformenergy},EquivGasSurfSpecsNames))...
            -sum(Matrix_H_ref_gas_species.*coeff_form_energy(:,indexformenergy)'))*0.0433634; % 1 kcal/mol = 0.0433634 eV
    end
end

for indexprintformenergy = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,['  ' ,sprintf('%1.3f', (FormEnergy(indexprintformenergy)))]);
end
fprintf (fileZacrosID_Simulation,[' # eV' CRLF]);


fprintf (fileZacrosID_Simulation,['gas_molec_weights     ']);
for indexweights = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,['  ' ,sprintf('%1.3f ', (MW(SpecIndx(GasSpecies(indexweights),EquivGasSurfSpecsNames)))*1000)]);
end
fprintf (fileZacrosID_Simulation,[' # g/mol' CRLF]);

fprintf (fileZacrosID_Simulation,['gas_molar_fracs         ']);
for indexpressures = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,['  ' ,sprintf('%1.3f ', (Molecular_Fractions(SpecIndx(GasSpecies(indexpressures),GasSpecies))))]);
end
fprintf (fileZacrosID_Simulation,[' # random value' CRLF]);
fprintf (fileZacrosID_Simulation,[CRLF]);

fprintf (fileZacrosID_Simulation,['n_surf_species             ', num2str(length(SurfaceSpecies)) CRLF]);

fprintf (fileZacrosID_Simulation,['surf_specs_names            ']);
for indexgasspecies = 1:length(SurfaceSpecies)
    fprintf (fileZacrosID_Simulation,[cell2mat(SurfaceSpecies(indexgasspecies)),'  ']);
end
fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['surf_specs_dent              ']);
for indexgasspecies = 1:length(SurfaceSpecies)
    fprintf (fileZacrosID_Simulation,[sprintf('%1.0f    ',1)]);
end
fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['snapshots                   on time  ', sprintf('%0.0e', 1e0) CRLF]);
fprintf (fileZacrosID_Simulation,['process_statistics          on time  ', sprintf('%0.0e', 1e0) CRLF]);
fprintf (fileZacrosID_Simulation,['species_numbers             on time  ', sprintf('%0.0e', 1e0) CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['# event_report                on' CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['max_steps                   infinity' CRLF]);
fprintf (fileZacrosID_Simulation,['max_time                    ' , sprintf('%1.1f', 200)  CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['wall_time                   ' , sprintf('%1.0f', 30), ' # in seconds'  CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,[' no_restart' CRLF]);

fprintf (fileZacrosID_Simulation,['' CRLF]);

fprintf (fileZacrosID_Simulation,['# debug_report_processes' CRLF]);
fprintf (fileZacrosID_Simulation,['# debug_report_global_energetics' CRLF]);
fprintf (fileZacrosID_Simulation,['# debug_check_processes' CRLF]);
fprintf (fileZacrosID_Simulation,['# debug_check_lattice' CRLF]);

fprintf (fileZacrosID_Simulation,[CRLF]);

fprintf (fileZacrosID_Simulation,['finish' CRLF]);

fclose(fileZacrosID_Simulation);

%From here we create the lattice_input file

fileZacrosID_Lattice = fopen('lattice_input.dat','w');

fprintf (fileZacrosID_Lattice,['lattice default_choice' CRLF]);

fprintf (fileZacrosID_Lattice,[CRLF]);

fprintf (fileZacrosID_Lattice,['hexagonal_periodic  ', sprintf('%1.1f %1.0f %1.0f ', 1 , 50 , 50), ...
    '# 1) lattice constant 2)copies cell horizontal direction 3)copies cell vertical direction' CRLF]);

fprintf (fileZacrosID_Lattice,[CRLF]);

fprintf (fileZacrosID_Lattice,['end_lattice' CRLF]);

fclose(fileZacrosID_Lattice);

%From here we create the energetics_input file

fileZacrosID_Energetics = fopen('energetics_input_from_Matlab.dat','w');

fprintf (fileZacrosID_Energetics,['energetics' CRLF]);

fprintf (fileZacrosID_Energetics,[CRLF]);
fprintf (fileZacrosID_Energetics,['#########################################################' CRLF]);

for indexenergetics = 1:length(GasSpecies)
    if indexenergetics <= length(Reference_Gas_Species)
        cluster_energy{indexenergetics} = - SurfSpecsQT(SpecIndx(Reference_Gas_Species{indexenergetics},Species_for_Q_for_energetics));
    else
        cluster_energy{indexenergetics} = FormEnergy(indexenergetics) - ...
            (SurfSpecsQT(SpecIndx(GasSpecies{indexenergetics},Species_for_Q_for_energetics)))*0.0433634; % 1 kcal/mol = 0.0433634 eV
    end
end

for indexcluster = 1:length(GasSpecies)
    fprintf (fileZacrosID_Energetics,[CRLF]);
    fprintf (fileZacrosID_Energetics,['cluster ', char(GasSpecies{indexcluster}), '_Point' CRLF]);
    fprintf (fileZacrosID_Energetics,['  sites 1 ' CRLF]);
    fprintf (fileZacrosID_Energetics,['  lattice_state ' CRLF]);
    fprintf (fileZacrosID_Energetics,['    1 ', char(GasSpecies{indexcluster}), '*   1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  site_types StTp1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  graph_multiplicity 1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  cluster_eng ' , sprintf('%1.3f', cluster_energy{indexcluster}),' # eV' CRLF]);
    fprintf (fileZacrosID_Energetics,['end_cluster' CRLF]);
    fprintf (fileZacrosID_Energetics,['' CRLF]);
    fprintf (fileZacrosID_Energetics,['#########################################################' CRLF]);
end

fprintf (fileZacrosID_Energetics,['' CRLF]);
fprintf (fileZacrosID_Energetics,['end_energetics' CRLF]);

fclose(fileZacrosID_Energetics);
