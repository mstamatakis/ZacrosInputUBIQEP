clear all
close all
clc

global Rgas GasFormH gammaRh BondIndex A_s_Forw A_s_Back T T0 ...
    beta_Forw beta_Back GasSpecies SurfaceSpecies MW QT

Rgas = 8.314/4.184/1000; %[kcal/mol/K]
conversion_kcal_over_mol_to_eV = 0.0433634; % 1 kcal/mol = 0.0433634 eV
Navo = 6.022e23; %[molec/mol]
h = 6.626070040e-34; %[Js]	

T0 = 300; %[K]
T = 500 + 273.15; %[K]
deltaT = T - T0; %[K]

gammaRh = 2.49e-8; %kmol/m2

GasSpecies = {'CH4','H2','H2O','CO','CO2'};
SurfaceSpecies = {'H*','OH*','H2O*','CO*','CO2*','C*','CH*','CH2*','CH3*','COOH*','O*'};

%Insert molecular weight of surface species
MW = zeros(1,(length(SurfaceSpecies)+length(GasSpecies)));
MW(SpecIndx('H*',SurfaceSpecies)) = 1.00794 /1000; %[kg/mol]
MW(SpecIndx('OH*',SurfaceSpecies)) = 17.00734 /1000; %[kg/mol]
MW(SpecIndx('H2O*',SurfaceSpecies)) = 18.01528 /1000; %[kg/mol]
MW(SpecIndx('CO*',SurfaceSpecies)) = 28.0101 /1000; %[kg/mol]
MW(SpecIndx('CO2*',SurfaceSpecies)) = 44.0095 /1000; %[kg/mol]
MW(SpecIndx('C*',SurfaceSpecies)) = 12.0107 /1000; %[kg/mol]
MW(SpecIndx('CH*',SurfaceSpecies)) = 13.01864 /1000; %[kg/mol]
MW(SpecIndx('CH2*',SurfaceSpecies)) = 14.02658 /1000; %[kg/mol]
MW(SpecIndx('CH3*',SurfaceSpecies)) = 15.9994 /1000; %[kg/mol]
MW(SpecIndx('COOH*',SurfaceSpecies)) = 45.01744 /1000; %[kg/mol]
MW(SpecIndx('O*',SurfaceSpecies)) = 16.04246 /1000; %[kg/mol]

%Insert molecular weight of gas species
MW(length(SurfaceSpecies)+SpecIndx('CH4',GasSpecies)) =  16.04246 /1000; %[kg/mol]
MW(length(SurfaceSpecies)+SpecIndx('H2',GasSpecies)) =  2.01588 /1000; %[kg/mol]
MW(length(SurfaceSpecies)+SpecIndx('H2O',GasSpecies)) =  18.01528 /1000; %[kg/mol]
MW(length(SurfaceSpecies)+SpecIndx('CO',GasSpecies)) =  28.0101 /1000; %[kg/mol]
MW(length(SurfaceSpecies)+SpecIndx('CO2',GasSpecies)) =  44.0095 /1000; %[kg/mol]

P = 0.1 ; %bar %random value
Molecular_Fractions = zeros(1,length(GasSpecies));
Molecular_Fractions(SpecIndx('CH4',GasSpecies)) = 0.01;
Molecular_Fractions(SpecIndx('H2',GasSpecies)) = 0;
Molecular_Fractions(SpecIndx('H2O',GasSpecies)) = 0.02;
Molecular_Fractions(SpecIndx('CO',GasSpecies)) = 0;
Molecular_Fractions(SpecIndx('CO2',GasSpecies)) = 0;
Partial_Pressure = P*Molecular_Fractions; %bar

%Let's calculate the formation energies of each species in gas phase
for i = 1:(length(SurfaceSpecies)+length(GasSpecies))
    if i <= length(SurfaceSpecies)
        GasFormH(i) = GasFormationEnthalpy(SurfaceSpecies{i}(1:end-1),T);
    else
        GasFormH(i) = GasFormationEnthalpy(GasSpecies{i-length(SurfaceSpecies)},T);
    end
end

%Let's create the QT0_ZeroCoverage vector
QT0_ZeroCoverage = zeros(1,(length(SurfaceSpecies)+length(GasSpecies)));
QT0_ZeroCoverage(SpecIndx('H*',SurfaceSpecies))= 62.3; %[kcal/mol]
QT0_ZeroCoverage(SpecIndx('OH*',SurfaceSpecies))= 70;
QT0_ZeroCoverage(SpecIndx('H2O*',SurfaceSpecies))= 10.8;
QT0_ZeroCoverage(SpecIndx('CO*',SurfaceSpecies))= 38.5;
QT0_ZeroCoverage(SpecIndx('CO2*',SurfaceSpecies))= 5.2;
QT0_ZeroCoverage(SpecIndx('C*',SurfaceSpecies))= 159;
QT0_ZeroCoverage(SpecIndx('CH*',SurfaceSpecies))= 151.2;
QT0_ZeroCoverage(SpecIndx('CH2*',SurfaceSpecies))= 109.3;
QT0_ZeroCoverage(SpecIndx('CH3*',SurfaceSpecies))= 42.4;
QT0_ZeroCoverage(SpecIndx('COOH*',SurfaceSpecies))= 62.2;
QT0_ZeroCoverage(SpecIndx('O*',SurfaceSpecies))= 100;
QT0_ZeroCoverage(length(SurfaceSpecies)+SpecIndx('CH4',GasSpecies))= 6;


%In order to create the matrix of IE add the dependency between the species for which you want calculate the Q (1 term)
%and the species that influenced as lateral interaction (teta) (2 term)
IE = zeros((length(SurfaceSpecies)+length(GasSpecies)),length(SurfaceSpecies));

IE(SpecIndx('H*',SurfaceSpecies),SpecIndx('H*',SurfaceSpecies)) = -2.5; 
IE(SpecIndx('H*',SurfaceSpecies),SpecIndx('CO*',SurfaceSpecies)) = -3.7;

IE(SpecIndx('OH*',SurfaceSpecies),SpecIndx('OH*',SurfaceSpecies)) = -25;
IE(SpecIndx('OH*',SurfaceSpecies),SpecIndx('H2O*',SurfaceSpecies)) = 25;
IE(SpecIndx('OH*',SurfaceSpecies),SpecIndx('O*',SurfaceSpecies)) = -33;

IE(SpecIndx('H2O*',SurfaceSpecies),SpecIndx('OH*',SurfaceSpecies)) = IE(SpecIndx('OH*',SurfaceSpecies),SpecIndx('H2O*',SurfaceSpecies));
IE(SpecIndx('H2O*',SurfaceSpecies),SpecIndx('H2O*',SurfaceSpecies)) = -4.5;

IE(SpecIndx('CO*',SurfaceSpecies),SpecIndx('H*',SurfaceSpecies)) = IE(SpecIndx('H*',SurfaceSpecies),SpecIndx('CO*',SurfaceSpecies));
IE(SpecIndx('CO*',SurfaceSpecies),SpecIndx('CO*',SurfaceSpecies)) = -15;

IE(SpecIndx('O*',SurfaceSpecies),SpecIndx('O*',SurfaceSpecies)) = -26;

%Let's create the f vector that contains the dependency on Temperature
f = zeros(1,(length(GasSpecies)+length(SurfaceSpecies)));
f(SpecIndx('H*',SurfaceSpecies))= 1.5;
f(SpecIndx('OH*',SurfaceSpecies))= 2;
f(SpecIndx('H2O*',SurfaceSpecies))= 2.5;
f(SpecIndx('CO*',SurfaceSpecies))= 2;
f(SpecIndx('CO2*',SurfaceSpecies))= 2;
f(SpecIndx('C*',SurfaceSpecies))= 1.5;
f(SpecIndx('CH*',SurfaceSpecies))= 2;
f(SpecIndx('CH2*',SurfaceSpecies))= 2.5;
f(SpecIndx('CH3*',SurfaceSpecies))= 2.5;
f(SpecIndx('COOH*',SurfaceSpecies))= 2.5;
f(SpecIndx('O*',SurfaceSpecies))= 1.5;
f(length(SurfaceSpecies)+SpecIndx('CH4',GasSpecies))= 2;

%Let's create the teta Matrix that contains the lateral interactions
teta(SpecIndx('H*',SurfaceSpecies))= 0;
teta(SpecIndx('OH*',SurfaceSpecies))= 0;
teta(SpecIndx('H2O*',SurfaceSpecies))= 0;
teta(SpecIndx('CO*',SurfaceSpecies))= 0;
teta(SpecIndx('CO2*',SurfaceSpecies))= 0;
teta(SpecIndx('C*',SurfaceSpecies))= 0;
teta(SpecIndx('CH*',SurfaceSpecies))= 0;
teta(SpecIndx('CH2*',SurfaceSpecies))= 0;
teta(SpecIndx('CH3*',SurfaceSpecies))= 0;
teta(SpecIndx('COOH*',SurfaceSpecies))= 0;
teta(SpecIndx('O*',SurfaceSpecies))= 0;

% Let's put all together and calculate the QT that is the Heat of Chemisorption of Surface Species at the temperature of the system
QT = QT0_ZeroCoverage + (IE * teta')' - f * Rgas * deltaT; %[kcal/mol]

% Insert here the reactions
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

%From this point on calculations and printing
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
    
    if UBIType(i)==1 || UBIType(i)==4
    fprintf (fileZacrosID_Mechanism,['  pre_expon ', num2str(PreExponForw{i}/(Rgas*4.184*T),'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  pe_ratio ', num2str(PreExponForw{i}/PreExponBack{i},'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  activ_eng ', num2str(EactForw{i}*conversion_kcal_over_mol_to_eV,'%1.3f')  CRLF]); %[eV]
    elseif UBIType(i)==2
    fprintf (fileZacrosID_Mechanism,['  pre_expon ', num2str(PreExponForw{i}/(Rgas*4.184*1000*1000*T),'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  pe_ratio ', num2str(PreExponForw{i}/PreExponBack{i},'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  activ_eng ', num2str(EactForw{i}*conversion_kcal_over_mol_to_eV,'%1.3f')  CRLF]); %[eV]    
    elseif UBIType(i)==5 || UBIType(i)==6
    fprintf (fileZacrosID_Mechanism,['  pre_expon ', num2str(PreExponForw{i}*gammaRh,'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  pe_ratio ', num2str(PreExponForw{i}/PreExponBack{i},'%1.2e')  CRLF]);
    fprintf (fileZacrosID_Mechanism,['  activ_eng ', num2str(EactForw{i}*conversion_kcal_over_mol_to_eV,'%1.3f')  CRLF]); %[eV]    
    end
    
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

%Diffusion
for indexdiff = 1:length(SurfaceSpecies)
    fprintf (fileZacrosID_Mechanism,['' CRLF]);
    fprintf (fileZacrosID_Mechanism,['reversible_step ', char(SurfaceSpecies{indexdiff}),'-diffusion' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  sites 2' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  neighboring 1-2' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  initial' CRLF]);
    fprintf (fileZacrosID_Mechanism,['    1 ',char(SurfaceSpecies{indexdiff}),'  1' CRLF]);
    fprintf (fileZacrosID_Mechanism,['    2 *  1' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  final' CRLF]);
    fprintf (fileZacrosID_Mechanism,['    1 *  1' CRLF]);
    fprintf (fileZacrosID_Mechanism,['    2 ',char(SurfaceSpecies{indexdiff}),'  1' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  site_types StTp1 StTp1' CRLF]);
    fprintf (fileZacrosID_Mechanism,['  pre_expon ',sprintf('%1.2e', (Rgas*4.184*1000)/Navo*T/h) CRLF]); %[1/s]
    fprintf (fileZacrosID_Mechanism,['  pe_ratio ',sprintf('%1.2e', 1) CRLF]); %for simmetry
    fprintf (fileZacrosID_Mechanism,['  activ_eng ',sprintf('%1.3f', 0.2*QT(SpecIndx(SurfaceSpecies{indexdiff},SurfaceSpecies))*conversion_kcal_over_mol_to_eV) CRLF]); 
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

random_seed = 71543 ;
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

%I chose as reference species CH4 H2 H2O

Reference_Gas_Species = {'CH4' , 'H2' , 'H2O'};

for indexRefGasSpecies = 1:length(Reference_Gas_Species)
    Matrix_H_ref_gas_species(indexRefGasSpecies) = GasFormH(length(SurfaceSpecies)+SpecIndx(Reference_Gas_Species{indexRefGasSpecies},GasSpecies));
end

MatrixCompRef = [1 0 0; 4 2 2; 0 0 1]; % in rows: C H O , in column: CH4, H2, H2O
vett_comp = zeros (length(Reference_Gas_Species),(length(SurfaceSpecies)+length(GasSpecies)));

%insert SurfaceSpaces here, if reference species put [0 0 0] , if not put components respect to [C H O]
vett_comp(:,SpecIndx('H*',SurfaceSpecies)) = [0 1 0]';
vett_comp(:,SpecIndx('OH*',SurfaceSpecies)) = [0 1 1]';
vett_comp(:,SpecIndx('H2O*',SurfaceSpecies)) = [0 2 1]';
vett_comp(:,SpecIndx('CO*',SurfaceSpecies)) = [1 0 1]';
vett_comp(:,SpecIndx('CO2*',SurfaceSpecies)) = [1 0 2]';
vett_comp(:,SpecIndx('C*',SurfaceSpecies)) = [1 0 0]';
vett_comp(:,SpecIndx('CH*',SurfaceSpecies)) = [1 1 0]';
vett_comp(:,SpecIndx('CH2*',SurfaceSpecies)) = [1 2 0]';
vett_comp(:,SpecIndx('CH3*',SurfaceSpecies)) = [1 3 0]';
vett_comp(:,SpecIndx('COOH*',SurfaceSpecies)) = [1 1 2]';
vett_comp(:,SpecIndx('O*',SurfaceSpecies)) = [0 0 1]';

%insert GasSpaces here, if reference species put [0 0 0] , if not put components respect to [C H O]
vett_comp(:,(length(SurfaceSpecies)+SpecIndx('CH4',GasSpecies))) = [0 0 0]';
vett_comp(:,(length(SurfaceSpecies)+SpecIndx('H2',GasSpecies)))= [0 0 0]';
vett_comp(:,(length(SurfaceSpecies)+SpecIndx('H2O',GasSpecies)))= [0 0 0]';
vett_comp(:,(length(SurfaceSpecies)+SpecIndx('CO',GasSpecies)))= [1 0 1]';
vett_comp(:,(length(SurfaceSpecies)+SpecIndx('CO2',GasSpecies)))= [1 0 2]';

for indexcomp = 1:length(vett_comp)
    coeff_form_energy(:,indexcomp) = inv(MatrixCompRef) * vett_comp(:,indexcomp);
end

for indexformenergy = 1:(length(SurfaceSpecies)+length(GasSpecies))
    if indexformenergy <= length(SurfaceSpecies)
        FormEnergy(indexformenergy) = (GasFormH(indexformenergy)-sum(Matrix_H_ref_gas_species.*coeff_form_energy(:,indexformenergy)'))*conversion_kcal_over_mol_to_eV; %[eV]
    else
        if sum(strcmp(GasSpecies{indexformenergy-length(SurfaceSpecies)},Reference_Gas_Species))==1
            FormEnergy(indexformenergy) = 0;
        elseif sum(strcmp(GasSpecies{indexformenergy-length(SurfaceSpecies)},Reference_Gas_Species))==0
            FormEnergy(indexformenergy) = (GasFormH(indexformenergy)-sum(Matrix_H_ref_gas_species.*coeff_form_energy(:,indexformenergy)'))*conversion_kcal_over_mol_to_eV; %[eV]
        end
    end
end

for indexprintformenergy = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,['  ' ,sprintf('%1.3f', (FormEnergy(length(SurfaceSpecies)+indexprintformenergy)))]);
end
fprintf (fileZacrosID_Simulation,[' # eV' CRLF]);


fprintf (fileZacrosID_Simulation,['gas_molec_weights     ']);
for indexweights = 1:length(GasSpecies)
    fprintf (fileZacrosID_Simulation,['  ' ,sprintf('%1.3f ', (MW(length(SurfaceSpecies)+SpecIndx(GasSpecies(indexweights),GasSpecies)))*1000)]);
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
coord_num = 6; %coordination number

fileZacrosID_Lattice = fopen('lattice_input.dat','w');

fprintf (fileZacrosID_Lattice,['lattice default_choice' CRLF]);

fprintf (fileZacrosID_Lattice,[CRLF]);

fprintf (fileZacrosID_Lattice,['hexagonal_periodic  ', sprintf('%1.1f %1.0f %1.0f ', 1 , 50 , 50), ...
    '# 1) lattice constant 2)copies cell horizontal direction 3)copies cell vertical direction' CRLF]);

fprintf (fileZacrosID_Lattice,[CRLF]);

fprintf (fileZacrosID_Lattice,['end_lattice' CRLF]);

fclose(fileZacrosID_Lattice);

%From here we create the energetics_input file
fileZacrosID_Energetics = fopen('energetics_input.dat','w');
fprintf (fileZacrosID_Energetics,['energetics' CRLF]);
fprintf (fileZacrosID_Energetics,[CRLF]);
fprintf (fileZacrosID_Energetics,['#########################################################' CRLF]);

for indexenergetics = 1:length(SurfaceSpecies)
    cluster_energy{indexenergetics} = FormEnergy(indexenergetics)-QT(indexenergetics)*conversion_kcal_over_mol_to_eV; %[eV]
end

for indexcluster = 1:length(SurfaceSpecies)
    fprintf (fileZacrosID_Energetics,[CRLF]);
    fprintf (fileZacrosID_Energetics,['cluster ', char(SurfaceSpecies{indexcluster}), '_Point' CRLF]);
    fprintf (fileZacrosID_Energetics,['  sites 1 ' CRLF]);
    fprintf (fileZacrosID_Energetics,['  lattice_state ' CRLF]);
    fprintf (fileZacrosID_Energetics,['    1 ', char(SurfaceSpecies{indexcluster}), '   1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  site_types StTp1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  graph_multiplicity 1' CRLF]);
    fprintf (fileZacrosID_Energetics,['  cluster_eng ' , sprintf('%1.3f', cluster_energy{indexcluster}),' # eV' CRLF]);
    fprintf (fileZacrosID_Energetics,['end_cluster' CRLF]);
    fprintf (fileZacrosID_Energetics,['' CRLF]);
    fprintf (fileZacrosID_Energetics,['#########################################################' CRLF]);
end

for indexpair_row = 1:length(SurfaceSpecies)
    for indexpair_column = indexpair_row:length(SurfaceSpecies)
        if IE(indexpair_row,indexpair_column)>0 || IE(indexpair_row,indexpair_column)<0
            fprintf (fileZacrosID_Energetics,[CRLF]);
            fprintf (fileZacrosID_Energetics,['cluster ', char(SurfaceSpecies{indexpair_row}),'/', char(SurfaceSpecies{indexpair_column}), '_pair' CRLF]);
            fprintf (fileZacrosID_Energetics,['  sites 2 ' CRLF]);
            fprintf (fileZacrosID_Energetics,['  neighboring 1-2' CRLF]);
            fprintf (fileZacrosID_Energetics,['  lattice_state ' CRLF]);
            fprintf (fileZacrosID_Energetics,['    1 ', char(SurfaceSpecies{indexpair_row}), ' 1'  CRLF]);
            fprintf (fileZacrosID_Energetics,['    2 ', char(SurfaceSpecies{indexpair_column}), ' 1'  CRLF]);
            fprintf (fileZacrosID_Energetics,['  site_types StTp1 StTp1' CRLF]);
            if strcmp(SurfaceSpecies{indexpair_row},SurfaceSpecies{indexpair_column})==1
                fprintf (fileZacrosID_Energetics,['  graph_multiplicity 2' CRLF]);
            else
                fprintf (fileZacrosID_Energetics,['  graph_multiplicity 1' CRLF]);
            end
            fprintf (fileZacrosID_Energetics,['  cluster_eng ' , sprintf('%1.3f', IE(indexpair_row,indexpair_column)/coord_num*conversion_kcal_over_mol_to_eV),' # eV' CRLF]);
            fprintf (fileZacrosID_Energetics,['end_cluster' CRLF]);
            fprintf (fileZacrosID_Energetics,['' CRLF]);
            fprintf (fileZacrosID_Energetics,['#########################################################' CRLF]);
        end
    end
end

fprintf (fileZacrosID_Energetics,['' CRLF]);
fprintf (fileZacrosID_Energetics,['end_energetics' CRLF]);

fclose(fileZacrosID_Energetics);
