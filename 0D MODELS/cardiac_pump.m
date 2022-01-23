% Model cardiac pump
% ------------------
%
% Order and names of the 24 variables in Vini
% 1 --> VLV: volume left ventricule
% 2 --> PLV: pression left ventricule
% 3 --> QMI: flowrate at the mitral valve
% 4 --> VLA: volume left atrium
% 5 --> PLA: pressure left atrium
% 6 --> QAO: flowrate at the aortic valve
% 7 --> PSAS: pression systemic aortic sinus
% 8 --> QSAS: flowrate systemic aortic sinus
% 9 --> PSAT: pressure systemic artery
% 10 --> QSAT: florate systemic artery
% 11 --> PSVN: pressure systemic vein
% 12 --> QSVN: flowrate systemic vein
%
% 13 --> VRV: volume right ventricule
% 14 --> PRV: pression right ventricule
% 15 --> QTI: flowrate at the tricuspid valve
% 16 --> VRA: volume right atrium
% 17 --> PRA: pressure right atrium
% 18 --> QPO: flowrate at the pulmonary aortic valve
% 19 --> PPAS: pression pulmonary artery sinus
% 20 --> QPAS: flowrate pulmonary artery sinus
% 21 --> PPAT: pressure pulmonary artery
% 22 --> QPAT: florate pulmonary artery
% 23 --> PPVN: pressure pulmonary vein
% 24 --> QPVN: flowrate pulmonary vein


            
function cardiac_pump()
        
    % Variables Initialisation and time period discretisation
    [Vini,tini,DT] = init();
    time = tini;
    
    % Paramètres de visualisation
    cvisu_ini=0;
    cvisu_fin=100;
    data1=1;
    data2=4;
    
    figure;
    for cycle =1:20
        [VSYS,VPUL] = solvePQV(Vini);
        if (cycle>cvisu_ini) & (cycle<cvisu_fin)
            plot(time,VSYS(data1,:),'k','LineWidth',3);
            hold on;
            plot(time,VSYS(data2,:),'r','LineWidth',3);
        end
        Vini = [VSYS(:,end)',VPUL(:,end)'];
        time = 1.0+time;
    end
    grid on;
    
    text_legend1 = leng_of_fig(data1);
    text_legend2 = leng_of_fig(data2);
    
    legend(text_legend1,text_legend2);
    hold off;
    

end

function [VARIABLE,t,DT] = init()
    
    % Time parameters
    [T,DT,Ts1,Ts2,Tpwb,Tpww] = time_parameters();
    % Heart parameters
    [LHEART_CONSTANT,RHEART_CONSTANT] = heart_data();
    % systemic & pulmonary parameters
    [SYS_CONSTANT,PUL_CONSTANT] = circulation_data();
    % Ventricular and atrium elastances
    [ELV,ERV,ELA,ERA,t] = elast(LHEART_CONSTANT,RHEART_CONSTANT,T,DT,Ts1,Ts2,Tpwb,Tpww);
    
    % Initialisation variables
    % Systemic circulation
    PSAS(1) = 80;
    VLV(1) = 120;
    VLA(1) = 60;
    QSAS(1) = 0;
    PSAT(1) = 80;
    QSAT(1) = 0;
    PSVN(1) = 10;
    % Pulmonary circulation
    PPAS(1) = 12;
    VRV(1) = 100;
    VRA(1) = 40;
    QPAS(1) = 0;
    PPAT(1) = 12;
    QPAT(1) = 0;
    PPVN(1) = 6;
        
    i = 1;
    % Systemic
    PLV(i) = LHEART_CONSTANT(5) + ELV(i)*(VLV(i)-LHEART_CONSTANT(6));
    PLA(i) =LHEART_CONSTANT(9) + ELA(i)*(VLA(i)-LHEART_CONSTANT(10));
    QAO(i) = VALVE_FLOW(PLV(i),PSAS(i),LHEART_CONSTANT(1));
    QMI(i) = VALVE_FLOW(PLA(i),PLV(i),LHEART_CONSTANT(2));
    % Pulmonary
    PRV(i) = RHEART_CONSTANT(5) + ERV(i)*(VRV(i)-RHEART_CONSTANT(6));
    PRA(i) =RHEART_CONSTANT(9) + ERA(i)*(VRA(i)-RHEART_CONSTANT(10));
    QPO(i) = VALVE_FLOW(PRV(i),PPAS(i),RHEART_CONSTANT(1));
    QTI(i) = VALVE_FLOW(PRA(i),PRV(i),RHEART_CONSTANT(2));
        
    QSVN(i) = (PSVN(i)-PRA(i))/SYS_CONSTANT(9);
    QPVN(i) = (PPVN(i)-PLA(i))/PUL_CONSTANT(9);

    VARIABLE = [VLV,PLV,QMI,VLA,PLA,QAO,PSAS,QSAS,PSAT,QSAT,PSVN,QSVN,...
                VRV,PRV,QTI,VRA,PRA,QPO,PPAS,QPAS,PPAT,QPAT,PPVN,QPVN];

        
end

function [VSYS,VPUL] = solvePQV(VINI)
    
    VLV(1) = VINI(1);
    PLV(1) = VINI(2);
    QMI(1) = VINI(3);
    VLA(1) = VINI(4);
    PLA(1) = VINI(5);
    QAO(1) = VINI(6);
    PSAS(1) = VINI(7);
    QSAS(1) = VINI(8);
    PSAT(1) = VINI(9);
    QSAT(1) = VINI(10);
    PSVN(1) = VINI(11);
    QSVN(1) = VINI(12);
    
    VRV(1) = VINI(1+12);
    PRV(1) = VINI(2+12);
    QTI(1) = VINI(3+12);
    VRA(1) = VINI(4+12);
    PRA(1) = VINI(5+12);
    QPO(1) = VINI(6+12);
    PPAS(1) = VINI(7+12);
    QPAS(1) = VINI(8+12);
    PPAT(1) = VINI(9+12);
    QPAT(1) = VINI(10+12);
    PPVN(1) = VINI(11+12);
    QPVN(1) = VINI(12+12);
    
    [T,DT,Ts1,Ts2,Tpwb,Tpww] = time_parameters();
    [SYS_CONSTANT,PUL_CONSTANT] = circulation_data();
    [LHEART_CONSTANT,RHEART_CONSTANT] = heart_data();
    [ELV,ERV,ELA,ERA,t] = elast(LHEART_CONSTANT,RHEART_CONSTANT,T,DT,Ts1,Ts2,Tpwb,Tpww);
    
    for i=1:length(t)-1
        
        VLV(i+1) = VLV(i) + DT*(QMI(i)-QAO(i));
        VLA(i+1) = VLA(i) + DT*(QPVN(i)-QMI(i));
        PSAS(i+1) = PSAS(i) + DT*(QAO(i)-QSAS(i))/SYS_CONSTANT(1);
        QSAS(i+1) = QSAS(i) + DT*(PSAS(i)-PSAT(i)-SYS_CONSTANT(2)*QSAS(i))/SYS_CONSTANT(3);
        PSAT(i+1) = PSAT(i) + DT*(QSAS(i)-QSAT(i))/SYS_CONSTANT(4);
        QSAT(i+1) = QSAT(i) + DT*(PSAS(i)-PSVN(i)-(SYS_CONSTANT(5)+SYS_CONSTANT(7)+SYS_CONSTANT(8))*QSAS(i))/SYS_CONSTANT(6);
        PSVN(i+1) = PSVN(i) + DT*(QSAT(i)-QSVN(i))/SYS_CONSTANT(10);
        
        VRV(i+1) = VRV(i) + DT*(QTI(i)-QPO(i));
        VRA(i+1) = VRA(i) + DT*(QSVN(i)-QTI(i));
        PPAS(i+1) = PPAS(i) + DT*(QPO(i)-QPAS(i))/PUL_CONSTANT(1);
        QPAS(i+1) = QPAS(i) + DT*(PPAS(i)-PPAT(i)-PUL_CONSTANT(2)*QPAS(i))/PUL_CONSTANT(3);
        PPAT(i+1) = PPAT(i) + DT*(QPAS(i)-QPAT(i))/PUL_CONSTANT(4);
        QPAT(i+1) = QPAT(i) + DT*(PPAS(i)-PPVN(i)-(PUL_CONSTANT(5)+PUL_CONSTANT(7)+PUL_CONSTANT(8))*QPAS(i))/PUL_CONSTANT(6);
        PPVN(i+1) = PPVN(i) + DT*(QPAT(i)-QPVN(i))/PUL_CONSTANT(10);

        % Systemic
        PLV(i+1) = LHEART_CONSTANT(5) + ELV(i+1)*(VLV(i+1)-LHEART_CONSTANT(6));
        PLA(i+1) =LHEART_CONSTANT(9) + ELA(i+1)*(VLA(i+1)-LHEART_CONSTANT(10));
        QAO(i+1) = VALVE_FLOW(PLV(i+1),PSAS(i+1),LHEART_CONSTANT(1));
        QMI(i+1) = VALVE_FLOW(PLA(i+1),PLV(i+1),LHEART_CONSTANT(2));
        % Pulmonary
        PRV(i+1) = RHEART_CONSTANT(5) + ERV(i+1)*(VRV(i+1)-RHEART_CONSTANT(6));
        PRA(i+1) = RHEART_CONSTANT(9) + ERA(i+1)*(VRA(i+1)-RHEART_CONSTANT(10));
        QPO(i+1) = VALVE_FLOW(PRV(i+1),PPAS(i+1),RHEART_CONSTANT(1));
        QTI(i+1) = VALVE_FLOW(PRA(i+1),PRV(i+1),RHEART_CONSTANT(2));
        
        QSVN(i+1) = (PSVN(i+1)-PRA(i+1))/SYS_CONSTANT(9);
        QPVN(i+1) = (PPVN(i+1)-PLA(i+1))/PUL_CONSTANT(9);
        
    end

    VSYS = [VLV;PLV;QMI;VLA;PLA;QAO;PSAS;QSAS;PSAT;QSAT;PSVN;QSVN];
    VPUL = [VRV;PRV;QTI;VRA;PRA;QPO;PPAS;QPAS;PPAT;QPAT;PPVN;QPVN];
end


function [ELV,ERV,ELA,ERA,t] = elast(LHEART_CONSTANT,RHEART_CONSTANT,T,DT,Ts1,Ts2,Tpwb,Tpww)

    % time
    t=[0:DT:T];
    
    % Variable Elastances
    ELV = LHEART_CONSTANT(4) + (LHEART_CONSTANT(3)-LHEART_CONSTANT(4))/2*vent_act(t,Ts1,Ts2);
    ERV = RHEART_CONSTANT(4) + (RHEART_CONSTANT(3)-RHEART_CONSTANT(4))/2*vent_act(t,Ts1,Ts2);
    ELA = LHEART_CONSTANT(8) + (LHEART_CONSTANT(7)-LHEART_CONSTANT(8))/2*atr_act(t,Tpwb,Tpww);
    ERA = RHEART_CONSTANT(8) + (RHEART_CONSTANT(7)-RHEART_CONSTANT(8))/2*atr_act(t,Tpwb,Tpww);

end


function [vact] = vent_act(t,Ts1,Ts2)
    % activation functions
    for i=1:length(t)
        if t(i)<Ts1
            vact(i) = 1-cos(t(i)*pi/Ts1);
        elseif t(i)<Ts2
            vact(i) = 1+cos( (t(i)-Ts1)*pi/(Ts2-Ts1) );
        else
            vact(i) = 0;
        end
    end
    
end

function [aact] = atr_act(t,Tpwb,Tpww)
    % activation functions
    for i=1:length(t)
        if t(i)<Tpwb
            aact(i) = 0;
        elseif t(i)<Tpwb+Tpww
            aact(i) = 1-cos( (t(i)-Tpwb)*2*pi/Tpww );
        else
            aact(i) = 0;
        end
    end
    
end

function [QV] = VALVE_FLOW(PL,PR,CQ)

    if PL<PR
        AR = 0.0;
    else
        AR = 1.0;
    end
    QV = CQ*AR*sqrt(abs(PL-PR));
    
end

function [LHEART_CONSTANT,RHEART_CONSTANT] = heart_data()
    
    %Left Heart constants
    CQao = 350;
    CQmi = 400;
    Elvs = 2.5;
    Elvd = 0.1;
    Plv0 = 1.0;
    Vlv0 = 5.0;
    Elamax = 0.25;
    Elamin = 0.15;
    Pla0 = 1.0;
    Vla0 = 4.0;
    LHEART_CONSTANT = [CQao CQmi Elvs Elvd Plv0 Vlv0 Elamax Elamin Pla0 Vla0];
    
    %Right Heart constants
    CQpo = 350;
    CQti = 400;
    Ervs = 1.15;
    Ervd = 0.1;
    Prv0 = 1.0;
    Vrv0 = 10.0;
    Eramax = 0.25;
    Eramin = 0.15;
    Pra0 = 1.0;
    Vra0 = 4.0;
    RHEART_CONSTANT = [CQpo CQti Ervs Ervd Prv0 Vrv0 Eramax Eramin Pra0 Vra0];
    
end

function [SYS_CONSTANT,PUL_CONSTANT] = circulation_data()
    
    % Systemic parameters
    Csas = 0.08;
    Rsas = 0.003;
    Lsas = 0.000062;
    Csat = 1.6;
    Rsat = 0.05;
    Lsat = 0.0017;
    Rsar = 0.5;
    Rscp = 0.52;
    Rsvn = 0.075;
    Csvn = 20.5;
    Csvc = 1.5;
    Vlv0 = 500;
    SYS_CONSTANT = [Csas Rsas Lsas Csat Rsat Lsat Rsar Rscp Rsvn Csvn Csvc Vlv0];
    
    % Pulmonary parameters
    Cpas = 0.18;
    Rpas = 0.002;
    Lpas = 0.000052;
    Cpat = 3.8;
    Rpat = 0.01;
    Lpat = 0.0017;
    Rpar = 0.05;
    Rpcp = 0.25;
    Rpvn = 0.006;
    Cpvn = 20.5;
    Cpvc = 1.5;
    Vrv0 = 400;
    PUL_CONSTANT = [Cpas Rpas Lpas Cpat Rpat Lpat Rpar Rpcp Rpvn Cpvn Cpvc Vrv0];
    
end

function [T,DT,Ts1,Ts2,Tpwb,Tpww] = time_parameters()

    T = 1;
    DT = 1.0E-04;
    Ts1 = 0.3;
    Ts2 = 0.45;
    Tpwb = 0.92;
    Tpww = 0.09;
    
end

function [ynp1]=rk4(yn,f,h)

    k1=f;
    k2=yn+h/2*k1;
    k3=yn+h/2*k2;
    k4=yn+h*k3;
    ynp1 = yn + h/6*(k1+2*k2+2*k3+k4);
    
end

function [text_legend] = leng_of_fig(data)

    switch data
        case 1
            text_legend = 'VLV: volume left ventricule';
        case 2
            text_legend = 'PLV: pression left ventricule';
        case 3
            text_legend = 'QMI: flowrate at the mitral valve';
        case 4
            text_legend = 'VLA: volume left atrium';
        case 5
            text_legend = 'PLA: pressure left atrium';
        case 6
            text_legend = 'QAO: flowrate at the aortic valve';
        case 7
            text_legend = 'PSAS: pression systemic aortic sinus';
        case 8
            text_legend = 'QSAS: flowrate systemic aortic sinus';
        case 9
            text_legend = 'PSAT: pressure systemic artery';
        case 10
            text_legend = 'QSAT: florate systemic artery';
        case 11
            text_legend = 'PSVN: pressure systemic vein';
        case 12
            text_legend = 'QSVN: flowrate systemic vein';
        case 13
            text_legend = 'VRV: volume right ventricule';
        case 14
            text_legend = 'PRV: pression right ventricule';
        case 15
            text_legend = 'QTI: flowrate at the tricuspid valve';
        case 16
            text_legend = 'VRA: volume right atrium';
        case 17
            text_legend = 'PRA: pressure right atrium';
        case 18
            text_legend = 'QPO: flowrate at the pulmonary aortic valve';
        case 19
            text_legend = 'PPAS: pression pulmonary artery sinus';
        case 20
            text_legend = 'QPAS: flowrate pulmonary artery sinus';
        case 21
            text_legend = 'PPAT: pressure pulmonary artery';
        case 22
            text_legend = 'QPAT: florate pulmonary artery';
        case 23
            text_legend = 'PPVN: pressure pulmonary vein';
        case 24
            text_legend = 'QPVN: flowrate pulmonary vein';
    end

end
