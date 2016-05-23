function out1 = EactUBIQEP(GasReact, SurfReact, GasProd, SurfProd, UBIType ,ReactionNumber)

global  SurfSpecsQT EquivGasSurfSpecsNames EqGasSurfaceSpecsFormH gammaRh BondIndex A_s_Forw A_s_Back T T0 beta_Forw beta_Back  Rgas MW

if UBIType == 1 % Non-activated atomic or non-dissociative molecular adsorption, e.g., A+*-->A*
    
    for i = 1:length(GasReact)
        HGasReact(i)=EqGasSurfaceSpecsFormH(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
        MWGasReact(i)=MW(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
    end
    
    for i = 1:length(SurfProd)
        HSurfProd(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
        QSurfProd(i)=SurfSpecsQT(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
    end
    
    Dab=0;
    deltaHrxn=-QSurfProd;
    EactForw = 0;
    out1(1) = max([0,deltaHrxn,EactForw]); %Eact forward effective
    out1(2) = out1(1) - deltaHrxn; %Eact backward effective
    n = 1; %order of reaction
    out1(3) = A_s_Forw(ReactionNumber)/gammaRh^n*sqrt((Rgas*1000*4.184)*T/(2*pi*MWGasReact))*(T/T0)^beta_Forw(ReactionNumber);
    out1(4) = A_s_Back(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Back(ReactionNumber);
    out1(5) = out1(3)*exp(-out1(1)/(Rgas*T));
    out1(6) = out1(4)*exp(-out1(2)/(Rgas*T));
    
end

if UBIType == 2 % Non-activated homonuclear dissociative adsorption, e.g., A2+2*-->2A*
    
    for i = 1:length(GasReact)
        HGasReact(i)=EqGasSurfaceSpecsFormH(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
        MWGasReact(i)=MW(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
    end
    
    for i = 1:length(SurfProd)
        HSurfProd(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
        QSurfProd(i)=SurfSpecsQT(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
    end
    
    Dab=sum(HSurfProd)-sum(HGasReact);
    deltaHrxn=Dab-sum(QSurfProd);
    EactForw = 0;
    out1(1) = max([0,deltaHrxn,EactForw]); %Eact forward effective
    out1(2) = out1 - deltaHrxn; %Eact backward effective
    n = 2; %order of reaction
    out1(3) = A_s_Forw(ReactionNumber)/gammaRh^n*sqrt((Rgas*1000*4.184)*T/(2*pi*MWGasReact))*(T/T0)^beta_Forw(ReactionNumber);
    out1(4) = A_s_Back(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Back(ReactionNumber);
    out1(5) = out1(3)*exp(-out1(1)/(Rgas*T));
    out1(6) = out1(4)*exp(-out1(2)/(Rgas*T));
end

if UBIType == 4 % Activated heteronuclear adsorption, e.g., AB+2*-->A*+B*
    
    for i = 1:length(GasReact)
        HGasReact(i)=EqGasSurfaceSpecsFormH(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
        QGasReact(i)=SurfSpecsQT(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
        MWGasReact(i)=MW(SpecIndx(GasReact{i},EquivGasSurfSpecsNames));
    end
    
    for i = 1:length(SurfProd)
        HSurfProd(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
        QSurfProd(i)=SurfSpecsQT(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
    end
    
    Dab=sum(HSurfProd)-sum(HGasReact);
    deltaHrxn=Dab-sum(QSurfProd);
    EactForw = BondIndex(ReactionNumber)*(deltaHrxn-QGasReact+prod(QSurfProd)/sum(QSurfProd));
    out1(1) = max([0,deltaHrxn,EactForw]); %Eact forward effective
    out1(2) = out1 - deltaHrxn; %Eact backward effective
    n = 2; %order of reaction
    out1(3) = A_s_Forw(ReactionNumber)/gammaRh^n*sqrt((Rgas*1000*4.184)*T/(2*pi*MWGasReact))*(T/T0)^beta_Forw(ReactionNumber);
    out1(4) = A_s_Back(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Back(ReactionNumber);
    out1(5) = out1(3)*exp(-out1(1)/(Rgas*T));
    out1(6) = out1(4)*exp(-out1(2)/(Rgas*T));
end

if UBIType == 5 % Heteronuclear surface dissociation, e.g., AB*+*-->A*+B*
    
    for i = 1:length(SurfReact)
        HSurfReact(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfReact{i},EquivGasSurfSpecsNames));
        QSurfReact(i)=SurfSpecsQT(SpecIndx(SurfReact{i},EquivGasSurfSpecsNames));
    end
    
    for i = 1:length(SurfProd)
        HSurfProd(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
        QSurfProd(i)=SurfSpecsQT(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
    end
    
    Dab=sum(HSurfProd)-sum(HSurfReact);
    deltaHrxn=Dab-sum(QSurfProd)+sum(QSurfReact);
    EactForw = BondIndex(ReactionNumber)*(deltaHrxn+prod(QSurfProd)/sum(QSurfProd));
    out1(1) = max([0,deltaHrxn,EactForw]); %Eact forward effective
    out1(2) = out1 - deltaHrxn; %Eact backward effective
    n = 2; %order of reaction
    out1(3) = A_s_Forw(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Forw(ReactionNumber);
    out1(4) = A_s_Back(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Back(ReactionNumber);
    out1(5) = out1(3)*exp(-out1(1)/(Rgas*T));
    out1(6) = out1(4)*exp(-out1(2)/(Rgas*T));
end

if UBIType == 6 % Surface disproportionation, e.g., A*+B*-->C*+D*
    
    for i = 1:length(SurfReact)
        HSurfReact(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfReact{i},EquivGasSurfSpecsNames));
        QSurfReact(i)=SurfSpecsQT(SpecIndx(SurfReact{i},EquivGasSurfSpecsNames));
    end
    
    for i = 1:length(SurfProd)
        HSurfProd(i)=EqGasSurfaceSpecsFormH(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
        QSurfProd(i)=SurfSpecsQT(SpecIndx(SurfProd{i},EquivGasSurfSpecsNames));
    end
    
    Dab=sum(HSurfProd)-sum(HSurfReact);
    deltaHrxn=Dab-sum(QSurfProd)+sum(QSurfReact);
    EactForw = BondIndex(ReactionNumber)*(deltaHrxn+prod(QSurfProd)/sum(QSurfProd));
    out1(1) = max([0,deltaHrxn,EactForw]); %Eact forward effective
    out1(2) = out1 - deltaHrxn; %Eact backward effective
    n = 2; %order of reaction
    out1(3) = A_s_Forw(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Forw(ReactionNumber);
    out1(4) = A_s_Back(ReactionNumber)/gammaRh^(n-1)*(T/T0)^beta_Back(ReactionNumber);
    out1(5) = out1(3)*exp(-out1(1)/(Rgas*T));
    out1(6) = out1(4)*exp(-out1(2)/(Rgas*T));
end


end
