% ------------------
% Model cardiac pump
% ------------------
% Korakianitis & Shi 2006
% NOT WORKING YET...
%
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


%% MAIN FUNCTION
            
function cardiac_pump_V2()

    global T Ts1 Ts2 Tpwb Tpww LHEART_CONSTANT RHEART_CONSTANT SYS_CONSTANT PUL_CONSTANT;
    
    % Time parameters
    [T,Ts1,Ts2,Tpwb,Tpww] = time_parameters();
    % Heart parameters
    [LHEART_CONSTANT,RHEART_CONSTANT] = heart_data();
    % systemic & pulmonary parameters
    [SYS_CONSTANT,PUL_CONSTANT] = circulation_data();

        
    % Variables Initialisation and time period discretisation
    Vini = init();
    
    % Paramètres de visualisation
    cvisu_ini=5;
    cvisu_fin=10;
    data1=1;
    data2=4;
    
    figure;
    
    for cycle =1:10
        [t,V] = solvePQV(Vini);
        time = t+(cycle-1)*T;
        if (cycle>=cvisu_ini) & (cycle<=cvisu_fin)
            plot(time,V(:,data1),'k','LineWidth',3);
            hold on;
            plot(time,V(:,data2),'r','LineWidth',3);
        end
        Vini(:) = V(end,:);
    end
    
    text_legend1 = leng_of_fig(data1);
    text_legend2 = leng_of_fig(data2);
    
    legend(text_legend1,text_legend2);
    hold off;
    

end

%% OTHER FUNCTIONS

function [VARIABLE] = init()
        
    % Initialisation variables
    % Systemic circulation
    PSAS = 100;
    VLV = 500;
    VLA = 20;
    QSAS = 0;
    PSAT = 100;
    QSAT = 0;
    PSVN = 0;
    % Pulmonary circulation
    PPAS =30;
    VRV = 500;
    VRA = 20;
    QPAS = 0;
    PPAT = 30;
    QPAT = 0;
    PPVN = 0;
    
    [PLV,PLA,QAO,QMI,PRV,PRA,QPO,QTI,QSVN,QPVN] = sys_pul(0,VLV,VLA,PSAS,PSVN,...
        VRV,VRA,PPAS,PPVN);
    
    VARIABLE = [VLV,PLV,QMI,VLA,PLA,QAO,PSAS,QSAS,PSAT,QSAT,PSVN,QSVN,...
                VRV,PRV,QTI,VRA,PRA,QPO,PPAS,QPAS,PPAT,QPAT,PPVN,QPVN];

        
end

function [PLV,PLA,QAO,QMI,PRV,PRA,QPO,QTI,QSVN,QPVN] = sys_pul(t,VLV,VLA,PSAS,PSVN,...
        VRV,VRA,PPAS,PPVN)
    
    global Ts1 Ts2 Tpwb Tpww LHEART_CONSTANT RHEART_CONSTANT SYS_CONSTANT PUL_CONSTANT;

    
    ELV = LHEART_CONSTANT(4) + (LHEART_CONSTANT(3)-LHEART_CONSTANT(4))/2*vent_act(t,Ts1,Ts2);
    ERV = RHEART_CONSTANT(4) + (RHEART_CONSTANT(3)-RHEART_CONSTANT(4))/2*vent_act(t,Ts1,Ts2);
    ELA = LHEART_CONSTANT(8) + (LHEART_CONSTANT(7)-LHEART_CONSTANT(8))/2*atr_act(t,Tpwb,Tpww);
    ERA = RHEART_CONSTANT(8) + (RHEART_CONSTANT(7)-RHEART_CONSTANT(8))/2*atr_act(t,Tpwb,Tpww);

    
    % Systemic
    PLV = LHEART_CONSTANT(5) + ELV*(VLV-LHEART_CONSTANT(6));
    PLA =LHEART_CONSTANT(9) + ELA*(VLA-LHEART_CONSTANT(10));
    QAO = VALVE_FLOW(PLV,PSAS,LHEART_CONSTANT(1));
    QMI = VALVE_FLOW(PLA,PLV,LHEART_CONSTANT(2));
    % Pulmonary
    PRV = RHEART_CONSTANT(5) + ERV*(VRV-RHEART_CONSTANT(6));
    PRA =RHEART_CONSTANT(9) + ERA*(VRA-RHEART_CONSTANT(10));
    QPO = VALVE_FLOW(PRV,PPAS,RHEART_CONSTANT(1));
    QTI = VALVE_FLOW(PRA,PRV,RHEART_CONSTANT(2));
        
    QSVN = (PSVN-PRA)/SYS_CONSTANT(9);
    QPVN = (PPVN-PLA)/PUL_CONSTANT(9);
        
end

function [t,V] = solvePQV(VINI)

    global T;
    
    VLV = VINI(1);
    PLV = VINI(2);
    QMI = VINI(3);
    VLA = VINI(4);
    PLA = VINI(5);
    QAO = VINI(6);
    PSAS = VINI(7);
    QSAS = VINI(8);
    PSAT = VINI(9);
    QSAT = VINI(10);
    PSVN = VINI(11);
    QSVN = VINI(12);
    
    VRV = VINI(1+12);
    PRV = VINI(2+12);
    QTI = VINI(3+12);
    VRA = VINI(4+12);
    PRA = VINI(5+12);
    QPO = VINI(6+12);
    PPAS = VINI(7+12);
    QPAS = VINI(8+12);
    PPAT = VINI(9+12);
    QPAT = VINI(10+12);
    PPVN = VINI(11+12);
    QPVN = VINI(12+12);
    
    Y0 = [VLV VLA PSAS QSAS PSAT QSAT PSVN VRV VRA PPAS QPAS PPAT QPAT PPVN];
    
    opts = odeset('RelTol',1e-06,'AbsTol',1e-06, 'MaxStep', 1);
    [t,Y] = ode15s(@fv,[0 T],Y0,opts);
    
    for i=1:length(t)
        VLV = Y(i,1);
        VLA = Y(i,2);
        PSAS = Y(i,3);
        PSVN = Y(i,7);
        VRV = Y(i,8);
        VRA = Y(i,9);
        PPAS = Y(i,10);
        PPVN = Y(i,14);
        [PLV,PLA,QAO,QMI,PRV,PRA,QPO,QTI,QSVN,QPVN] = sys_pul(t(i),VLV,VLA,PSAS,PSVN,...
        VRV,VRA,PPAS,PPVN);
        V(i,2) = PLV;
        V(i,3) = QMI;
        V(i,5) = PLA;
        V(i,6) = QAO;
        V(i,12) = QSVN;
        
        V(i,2+12) = PRV;
        V(i,3+12) = QTI;
        V(i,5+12) = PRA;
        V(i,6+12) = QPO;
        V(i,12+12) = QPVN;
    end

    V(:,1) = Y(:,1);
    V(:,4) = Y(:,2);
    V(:,7) = Y(:,3);
    V(:,8) = Y(:,4);
    V(:,9) = Y(:,5);
    V(:,10) = Y(:,6);
    V(:,11) = Y(:,7);

    V(:,1+12) = Y(:,8);
    V(:,4+12) = Y(:,9);
    V(:,7+12) = Y(:,10);
    V(:,8+12) = Y(:,11);
    V(:,9+12) = Y(:,12);
    V(:,10+12) = Y(:,13);
    V(:,11+12) = Y(:,14);
            
end

function F = fv(t,y)

    global SYS_CONSTANT PUL_CONSTANT;

    
   [~,~,QAO,QMI,~,~,QPO,QTI,QSVN,QPVN] = sys_pul(t,y(1),y(2),y(3),y(7),...
        y(8),y(9),y(10),y(14));
    
    F(1) = QMI-QAO;
    F(2) = QPVN-QMI;
    F(3) = (QAO-y(4))/SYS_CONSTANT(1);
    F(4) = (y(3)-y(5)-SYS_CONSTANT(2)*y(4))/SYS_CONSTANT(3);
    F(5) = (y(4)-y(6))/SYS_CONSTANT(4);
    F(6) = (y(3)-y(7)-(SYS_CONSTANT(5)+SYS_CONSTANT(7)+SYS_CONSTANT(8))*y(4))/SYS_CONSTANT(6);
    F(7) = (y(6)-QSVN)/SYS_CONSTANT(10);
    
    F(8) = QTI-QPO;
    F(9) = QSVN-QTI;
    F(10) = (QPO-y(11))/PUL_CONSTANT(1);
    F(11) = (y(10)-y(12)-PUL_CONSTANT(2)*y(11))/PUL_CONSTANT(3);
    F(12) = (y(11)-y(13))/PUL_CONSTANT(4);
    F(13) = (y(10)-y(14)-(PUL_CONSTANT(5)+PUL_CONSTANT(7)+PUL_CONSTANT(8))*y(11))/PUL_CONSTANT(6);
    F(14) = (y(13)-QPVN)/PUL_CONSTANT(10);
    
    F = F';
    
end


function [vact] = vent_act(t,Ts1,Ts2)
    % activation functions
    if t<Ts1
        vact = 1-cos(t*pi/Ts1);
    elseif t<Ts2
        vact = 1+cos( (t-Ts1)*pi/(Ts2-Ts1) );
    else
        vact = 0;
    end
    
end

function [aact] = atr_act(t,Tpwb,Tpww)
    % activation functions
    if t<Tpwb
        aact = 0;
    elseif t<Tpwb+Tpww
        aact = 1-cos( (t-Tpwb)*2*pi/Tpww );
    else
        aact = 0;
    end
    
end

function [QV] = VALVE_FLOW(PL,PR,CQ)

    if PL<=PR
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

function [T,Ts1,Ts2,Tpwb,Tpww] = time_parameters()

    T = 1;
    Ts1 = 0.3;
    Ts2 = 0.45;
    Tpwb = 0.92;
    Tpww = 0.09;
    
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
