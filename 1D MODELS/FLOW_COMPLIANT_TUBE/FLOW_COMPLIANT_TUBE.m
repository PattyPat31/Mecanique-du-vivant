%**************************************************
%
%  ECOULEMENT UNIDIMENSIONNEL DANS UN TUBE SOUPLE
%
%**************************************************
%
%
%              conditions aux limites 0D
%       ---------------------------------------
%
% ATTENTION: unités grammes, centimètres, secondes
%% MAIN FUNCTION

function FLOW_COMPLIANT_TUBE()

    global mmHg
    mmHg=1334;    % coefficient pour changer d'unités de pression dyn/cm2<-->mmHg
    
    % Fichier de sortie des données
    output_data = fopen('sortie.dat','w');

    %.........................................
    %        DONNEES DU TUBE
    %.......................................
    %
    FNAME = 'param.dat';
    [EL,A0,LT,RO,mu,ResT,Comp,DT,type,typs,Tcycle,LC,QMEAN,NMODE,frequencies,...
        AMP_flowrate,PHASE_flowrate,Pout,Pin,AMP_pressure,PHASE_pressure] = datatube(FNAME);
    c0=sqrt(EL/RO);   % vitesse propagation onde
    Z0=RO*c0/A0;  % impédance du tube
    % Décomposition de la résistance périphérique en 2 parties: Z0 et Res
    Res=ResT-Z0;
    
    % Paramètres numériques
    %***
    % calcul du nombre de points en x pour respecter le CFL>c0
    CFL = 1.1*c0;
    PAS_X_MIN = CFL*DT;
    N_MAX = floor(LT/PAS_X_MIN) + 1;
    % Nombre de points en x
    if mod(N_MAX,2)==0 
      NX = N_MAX+1; % cas nombre pair
    else
      NX = N_MAX; % cas nombre impair
    end
    DX=LT/(NX-1);   % Pas en x
    x = [0:DX:LT];
    NTIME=Tcycle/DT; % Nombre de points en temps dans un cycle
    R=DT/DX;  % Inverse du CFL
    R2=R/2;
    
    % Autres paramètres numériques
    %***
    NBB=1;
    I1=1;
    I2=2;
    I3=NX-1;
    I4=NX;
    RPI=sqrt(pi);

    % INITIALISATION pression
%     PEM= Pin;
%     for mode=1:NMODE
%         PEM = PEM + AMP_pressure(mode)*mmHg*cos(PHASE_pressure(mode));
%     end
%     Pin/mmHg
%     Pout/mmHg
%     PEM/mmHg
%     (ResT*QMEAN+Pout)/mmHg
%     pause
    PEM=Pin;
    PSM=PEM-8*pi*mu*NBB/RO*QMEAN*LT/A0^2;
    DPI=(PEM-PSM)/(NX-1);
    PC=PSM-Z0*QMEAN;
    P(1)=PEM;
    P(NX)=PSM;
    for I=2:NX-1
       P(I) = P(I-1)-DPI;
    end

    % INITIALISATION section, vitesse, contrainte pariétale
    for I=1:NX
       A(I)=A0*(P(I)/EL+1);   % Loi du tube
       U(I)=QMEAN/A(I);   % Débit moyen
       T(I)=-4*mu*sqrt(pi/A(I))*U(I); % approximation Poiseuille de la contrainte pariétale.
    end

    ITC=0.0;

    %          CALCUL LC CYCLE(S)
    time(1) = 0.0;
    tsortie = 0.0;
    % -------------------
    % Début boucle cycles
    % -------------------
    for NC=1:LC

        dis = sprintf('Cycle N° %i ',NC);
        disp(dis);

        % Sauvegarde Data initialisation
        if (ITC==0)
           for I=1:NX
              DEB=A(I)*U(I);
              fprintf(output_data,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',DEB, A(I), U(I),P(I),T(I));
            end
        end

        ITC = ITC+1;

        % Paramètres pour les sorties data   
        MC=1;
        DMC=1;

        % Début boucle en temps (periode de Tcycle seconde)
        for N=1:NTIME

            time(NTIME*(NC-1)+N) = (NC-1)*Tcycle+N*DT;


            %.....................................
            %          CONDITIONS AUX LIMITES
            %.....................................
            % ENTREE
            switch type
                case 3
                % Pression à l'entrée du tube (sinusoidale)
                PLE=PEM+20*mmHg*(sin(2*pi*N*DT/Tcycle));
                ALE=(PLE/EL+1)*A0;
                S1=U(1)*U(1)/2.+P(1)/RO;
                S2=U(2)*U(2)/2.+P(2)/RO;
                TO=(2*RPI*DT/RO)*T(1)/sqrt(A(1));
                ULE=U(1)+R*(S1-S2)+TO;
                TLE = -4*mu*sqrt(pi/ALE)*ULE;

                case 4
                % Pression à l'entrée du tube (physiologique)
                PLE= Pin;
                for mode=1:NMODE
                    PLE = PLE + AMP_pressure(mode)*mmHg*cos(2*pi*frequencies(mode)*N*DT+PHASE_pressure(mode));
                end
                ALE=(PLE/EL+1)*A0;
                S1=U(1)*U(1)/2.+P(1)/RO;
                S2=U(2)*U(2)/2.+P(2)/RO;
                TO=(2*RPI*DT/RO)*T(1)/sqrt(A(1));
                ULE=U(1)+R*(S1-S2)+TO;
                TLE = -4*mu*sqrt(pi/ALE)*ULE;

                case 1
                % Débit à l'entrée du tube (sinusoidale)
                QLE=QMEAN+2*(sin(2*pi*N*DT/Tcycle));
                QLE=real(QLE);
                S1=U(1)*A(1);
                S2=U(2)*A(2);
                ALE=A(1)+R*(S1-S2);
                ULE=QLE/ALE;
                PLE=EL*(ALE/A0-1);
                TLE = -4*mu*sqrt(pi/ALE)*ULE;

                case 2
                % Débit à l'entrée du tube (physiologique)
                QLE= QMEAN;
                for mode=1:NMODE
                    QLE = QLE + real(AMP_flowrate(mode)*exp(1i*2*pi*frequencies(mode)*N*DT+1i*PHASE_flowrate(mode)));
                end
                QLE=real(QLE);
                S1=U(1)*A(1);
                S2=U(2)*A(2);
                ALE=A(1)+R*(S1-S2);
                ULE=QLE/ALE;
                PLE=EL*(ALE/A0-1);
                TLE = -4*mu*sqrt(pi/ALE)*ULE;                
            end
            
            % SORTIE
            switch typs
                case 2
                % Condition de sortie Windkessel de type RCR non réflexive
                DEB=A(NX)*U(NX);
                DEB1=A(NX-1)*U(NX-1);
                ALS=A(NX)-DT/DX*(DEB-DEB1);   % A*(n+1)
                PLS=EL*(ALS/A0-1);            % P*(n+1)
                QLSe=(PLS-PC)/Z0;      % Q*(n) débit dans R1
                QLS=(PC-Pout)/Res;        % Q(n) débit dans R2
                PC = PC+DT/Comp*(QLSe-QLS); % PC(n+1)
                QLSe=(PLS-PC)/Z0;        % Q*(n+1)
                QLS=(PC-Pout)/Res;        % Q(n+1)
                ULS=QLSe/ALS;
                TLS = -4*mu*sqrt(pi/ALS)*ULS;

                case 1
                % Condition de sortie Windkessel simple R
                DEB=A(NX)*U(NX);
                DEB1=A(NX-1)*U(NX-1);
                ALS=A(NX)-DT/DX*(DEB-DEB1);
                PLS=EL*(ALS/A0-1);             
                QLS=(PLS-0.0)/ResT;        % Q(n+1) débit dans R
                ULS=QLS/ALS;
                TLS = -4*mu*sqrt(pi/ALS)*ULS;
                
                case 3
                % Condition débit imposé
                QLS = QMEAN;
                S1=U(NX-1)*A(NX-1);
                S2=U(NX)*A(NX);
                ALS=A(NX)+R*(S1-S2);
                ULS=QLS/ALS;
                PLS=EL*(ALS/A0-1);
                TLS = -4*mu*sqrt(pi/ALS)*ULS;

                case 4
                % Condition pression imposée
                PLS=PSM;
                ALS=(PLS/EL+1)*A0;
                S1=U(NX-1)*U(NX-1)/2.+P(NX-1)/RO;
                S2=U(NX)*U(NX)/2.+P(NX)/RO;
                TO=(2*RPI*DT/RO)*T(NX)/sqrt(A(NX));
                ULS=U(NX)+R*(S1-S2)+TO;
                TLS = -4*mu*sqrt(pi/ALS)*ULS;
                
            end

            %......................................
            %
            %         CALCUL POINTS INTERNES
            %
            %.......................................

            %        CALCUL 1ER DEMI-PAS
            for I=I1:I3
              AA=A(I);
              UU=U(I);
              PP=P(I);
              TT=T(I);
              AA1=A(I+1);
              UU1=U(I+1);
              PP1=P(I+1);
              TT1=T(I+1);

              AM=0.5*(AA+AA1);
              UM=0.5*(UU+UU1);
              TM=0.5*(TT+TT1);

              A1Z=R2*(AA1*UU1-AA*UU);
              A1(I+1)=(AA1+AA)/2.-A1Z;
              Q1=UU1*UU1/2.+PP1/RO;
              Q=UU*UU/2.+PP/RO;
              Q2=(RPI*DT/RO)*NBB*(TM/sqrt(AM));
              U1(I+1)=(UU1+UU)/2.-R2*(Q1-Q)+Q2;
              P1(I+1)=EL*(A1(I+1)/A0-1);
              T1(I+1) = -4*mu*sqrt(pi/A1(I+1))*U1(I+1);
            end

            %        CALCUL 2EME DEMI-PAS
            for I=I2:I3
              AA=A1(I);
              UU=U1(I);
              PP=P1(I);
              TT=T1(I);
              AA1=A1(I+1);
              UU1=U1(I+1);
              PP1=P1(I+1);
              TT1=T1(I+1);

              AM=0.5*(AA+AA1);
              UM=0.5*(UU+UU1);
              TM=0.5*(TT+TT1);

              UN=U(I);

              AZ=-R*(AA1*UU1-AA*UU);
              A(I)=A(I)+AZ;
              S1=UU1*UU1/2.+PP1/RO;
              S=UU*UU/2.+PP/RO;
              UZ=-R*(S1-S)+(2*RPI*DT/RO)*(TM/sqrt(AM))*NBB;
              U(I)=U(I)+UZ;
              P(I)=EL*(A(I)/A0-1);
              T(I) = -4*mu*sqrt(pi/A(I))*U(I);
            end

            A(1)=ALE;
            U(1)=ULE;
            P(1)=PLE;
            T(1)=TLE;
            A(NX)=ALS;
            U(NX)=ULS;
            P(NX)=PLS;
            T(NX)=TLS;

          % Sortie data du 1er au LCème cycle tous les DMC pas de temps
          if ((ITC>=1)&&(MC==N))
            %fprintf(output_data,'%4i\n',N);
            tsortie = [tsortie time(NTIME*(NC-1)+N)];
            MC=MC+DMC;
            for I=1:NX
                DEB=A(I)*U(I);
                fprintf(output_data,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',DEB, A(I), U(I), P(I), T(I));
            end
          end

        end
        % Fin boucle en temps (periode de Tcycle secondes)


    end
    % -------------------
    % Fin boucle cycles
    % -------------------

    fclose(output_data);
    visu(DX,LT,NX,tsortie);

end
%% FUNCTIONS

% ----------------------------------------
% LECTURE DES DONNEES ET CALCUL PARAMETRES
% ----------------------------------------
function [EL,A0,LT,RO,mu,ResT,Comp,DT,type,typs,Tcycle,LC,...
    QMEAN,NMODE,frequencies,AMP_flowrate,PHASE_flowrate,Pout,Pin,...
    AMP_pressure,PHASE_pressure] = datatube(FNAME)

    global mmHg
    
    datafile = fopen(FNAME,'r');
    % Paramètres tube -- loi du tube: A = A0(P/EL+1)
    %***
    fscanf(datafile,'%s',4);
    EL = fscanf(datafile,'%f',1);  % Elastance
    fscanf(datafile,'%s',4);
    A0 = fscanf(datafile,'%f',1);  % Section (cm^2)
    fscanf(datafile,'%s',4);
    LT = fscanf(datafile,'%f',1);  % longueur du tube (cm)
    fscanf(datafile,'%s',4);
    RO = fscanf(datafile,'%f',1);  % densité du fluide
    fscanf(datafile,'%s',4);
    mu = fscanf(datafile,'%f',1);  % viscosité du fluide
    fscanf(datafile,'%s',3);
    ResT = fscanf(datafile,'%f',1);% Résistance périphérique totale
    fscanf(datafile,'%s',3);
    Comp = fscanf(datafile,'%f',1);% Compliance périphérique
    fscanf(datafile,'%s',4);
    DT = fscanf(datafile,'%f',1);  % Pas de temps
    fscanf(datafile,'%s',5);
    % Type de signal d'entrée
    %1: débit sinusoidal; 2: débit physiologique; 3: pression sinusoidale
    type = fscanf(datafile,'%i',1);  
    fscanf(datafile,'%s',5);
    % Type de condition de sortie
    %1: 0D Résistif; 2: 0D non reflexif RCR
    typs = fscanf(datafile,'%i',1);  
    fscanf(datafile,'%s',2);
    Tcycle = fscanf(datafile,'%f',1);% Pas de temps
    fscanf(datafile,'%s',4);
    LC = fscanf(datafile,'%i',1);    % nombre de cycle
    fscanf(datafile,'%s',6);
    % Débit moyen
    QMEAN = fscanf(datafile,'%f',1);
    % Nombre de modes de Fourier (pour la reconstruction des conditions aux limites
    % de débit/pression physiologique)
    NMODE = fscanf(datafile,'%i',1);  
    % Fréquences des modes
    frequencies = fscanf(datafile,'%f',NMODE);
    % Amplitudes des modes de débit
    AMP_flowrate = fscanf(datafile,'%f',NMODE);
    % Phase des modes de débit
    PHASE_flowrate = fscanf(datafile,'%f',NMODE);
    % Pression de sortie périphérique
    Pout = fscanf(datafile,'%f',1);
    Pout = Pout*mmHg;
    fscanf(datafile,'%s',3);
    % Pression moyenne d'entrée
    Pin = fscanf(datafile,'%f',1);
    Pin = Pin*mmHg;
    % Amplitudes des modes de pression
    AMP_pressure = fscanf(datafile,'%f',NMODE);
    % Phase des modes de pression
    PHASE_pressure = fscanf(datafile,'%f',NMODE);
end

% -----------------------
% VISUALISATION RESULTATS
% -----------------------
function []=visu(DX,LT,NX,tsortie)

    global mmHg

    x = [0:DX:LT];
    % Définir si nécessaire un indice supplémentaire de visualisation dans
    % le tube, comme par exemple: milieu = fix((NX-1)/2);
    
    output = fopen('sortie.dat','r');

    % Evolution en temps à l'entrée, au milieu et à la sortie du tube
    for i=1:length(tsortie)
        scan = fscanf(output,'%f',5);
        debit_entree(i) = scan(1);
        section_entree(i)= scan(2);
        velocity_entree(i) = scan(3);
        pressure_entree(i) = scan(4)/mmHg;
        tau_entree(i) = scan(5);

%         fscanf(output,'%f',5*(milieu-2));
%         scan = fscanf(output,'%f',5);
%         debit_milieu(i) = scan(1);
%         section_milieu(i)= scan(2);
%         velocity_milieu(i) = scan(3);
%         pressure_milieu(i) = scan(4)/mmHg;
%         tau_milieu(i) = scan(5);
% 
%         fscanf(output,'%f',5*(NX-milieu-1));

        fscanf(output,'%f',5*(NX-2));
        scan = fscanf(output,'%f',5);
        debit_sortie(i) = scan(1);
        section_sortie(i)= scan(2);
        velocity_sortie(i) = scan(3);
        pressure_sortie(i) = scan(4)/mmHg;
        tau_sortie(i) = scan(5);
    end

    % Enlever les commentaires suivant les figures souhaitées
    %--------------------------------------------------------
    %title1 = strcat('\bf Section time evolution (Début du tube)');
    %figure;plot(tsortie,section_entree);title(title1);xlabel('\bf time');grid on;
    %title2 = strcat('\bf Velocity time evolution (Début du tube)');
    %figure;plot(tsortie,velocity_entree);title(title2);xlabel('\bf time');grid on;
    title3 = strcat('\bf Pressure time evolution (Début du tube)');
    figure;plot(tsortie,pressure_entree);title(title3);xlabel('\bf time');grid on;
    %title4 = strcat('\bf Tau time evolution (Début du tube)');
    %figure;plot(tsortie,tau_entree);title(title4);xlabel('\bf time');grid on;
    title5 = strcat('\bf Debit time evolution (Début du tube)');        
    figure;plot(tsortie,debit_entree);title(title5);xlabel('\bf time');grid on;
    %close all;

    %title1 = strcat('\bf Section time evolution (Fin du tube)');
    %figure;plot(tsortie,section_sortie);title(title1);xlabel('\bf time');grid on;
    %title2 = strcat('\bf Velocity time evolution (Fin du tube)');
    %figure;plot(tsortie,velocity_sortie);title(title2);xlabel('\bf time');grid on;
    title3 = strcat('\bf Pressure time evolution (Fin du tube)');
    figure;plot(tsortie,pressure_sortie);title(title3);xlabel('\bf time');grid on;
    %title4 = strcat('\bf Tau time evolution (Fin du tube)');
    %figure;plot(tsortie,tau_sortie);title(title4);xlabel('\bf time');grid on;
    title5 = strcat('\bf Debit time evolution (Fin du tube)');        
    figure;plot(tsortie,debit_sortie);title(title5);xlabel('\bf time');grid on;
    %close all;

    %title1 = strcat('\bf Section time evolution (Milieu du tube)');
    %figure;plot(tsortie,section_milieu);title(title1);xlabel('\bf time');grid on;
    %title2 = strcat('\bf Velocity time evolution (Milieu du tube)');
    %figure;plot(tsortie,velocity_milieu);title(title2);xlabel('\bf time'); grid on;%axis([0 4 min(min(velocity_milieu)) max(max(velocity_milieu))]);grid on;
    %title3 = strcat('\bf Pressure time evolution (Milieu du tube)');
    %figure;plot(tsortie,pressure_milieu);title(title3);xlabel('\bf time'); grid on;%axis([0 4 min(min(pressure_milieu)) max(max(pressure_milieu))]);grid on;
    %title4 = strcat('\bf Tau time evolution (Milieu du tube)');
    %figure;plot(tsortie,tau_milieu);title(title4);xlabel('\bf time');grid on;
    %title5 = strcat('\bf Debit time evolution (Milieu du tube)');        
    %figure;plot(tsortie,debit_milieu);title(title5);xlabel('\bf time');grid on;

    fclose(output);

end




