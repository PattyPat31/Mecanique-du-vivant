% -------------------------------------------------------------------------
%
%               -------------------------------
%               MODEL 1D OF NETWORK BLOOD FLOW
%               -------------------------------
%
%                       NETWORKFLOW
%                           ***
%                      P. Cathalifaud
%
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%                       MAIN FUNCTION
% -------------------------------------------------------------------------

function NETWORKFLOW()

    % --------------------
    %     Parameters
    % --------------------
    
    % Fluid parameters
    RO=1.06;     
    ETA=0.035;

    % Network data reading
    FNAME = 'donnees.dat';    
    datafile = fopen(FNAME,'r');
    [NBR,DAD,SECTI,ELASTA,LONG,NB,NBPTS,entree,sortie,T_period,PSORTIE,NCycle,DT1,Msave] = data(datafile);
    [PSING,nbif,njonc,nelarg,fuite,special,BW] = connectivities(datafile);
    fclose(datafile);

    % Entry pressure data
    [PRE,NPRE] = entry();
        
    % Numerical parameters
    DT2=DT1;
    e=ceil(T_period/DT1)-fix(T_period/DT1);
    C1=RO/(48*ETA*pi);
    C2=4*ETA*sqrt(pi);
    C3=RO/(4*sqrt(pi));
    M2=DT2*sqrt(pi)/RO;
    PAS = LONG./NBPTS;
    RAP = DT2./PAS;
    RAPI = PAS/DT2;
            
    % --------------------
    % Data Initialisations
    % --------------------
    [A,U,P,T] = init(entree,sortie,DAD,ELASTA,SECTI,LONG,NB,NBPTS,ETA,NBR,8*pi*ETA/RO);
    PE1=0;
    PIE=0;

    %               Save t=0 into output files
    % ---------------------------------------------------------------------
    output_inlet = fopen('A_U_P_DEB_T_inlet.dat','w');
    for br=1:NBR
        i = 1;
        DEB = A(br,i)*U(br,i);
        fprintf(output_inlet,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
    end
    output_pt4 = fopen('A_U_P_DEB_T_pt4.dat','w');
    for br=1:NBR
        i = 4;
        DEB = A(br,i)*U(br,i);
        fprintf(output_pt4,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
    end
    output_outlet = fopen('A_U_P_DEB_T_outlet.dat','w');
    for br=1:NBR
        i = NBPTS(br)+1;
        DEB = A(br,i)*U(br,i);
        fprintf(output_outlet,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
    end
    % ---------------------------------------------------------------------


    VOLT=0;
    for br=1:NBR
        IL=NBPTS(br)+1;
        VOLTI=0.5*LONG(br)*(A(br,1)+A(br,IL));
        VOLT=VOLT+VOLTI;
        dis = sprintf('branche N° %i; volume=%12.5f',br,VOLTI);
        disp(dis);
    end       
    dis = sprintf('Volume total = %12.5f',VOLT);
    disp(dis);

    
    
%   ***********************************
%       COMPUTATION OF THE FLOW
%   ***********************************

    for NC=1:NCycle  % number of cycles
        VVT=0;
        M=Msave;  % saving iteration

        IL=NBPTS(NBR)+1;
        QEI=sum(A(entree',1).*U(entree',1));
        QE0=QEI;
        QSI=sum(A(sortie',IL).*U(sortie',IL));
        if fuite>0
            QSI=QSI+A(fuite,1)*U(fuite,1);
        end
        QS0=QSI;
        
        % -------------------
        % BEGINNING TIME LOOP
        % -------------------
        for N=1:fix(T_period/DT1+e)
            AZT=0;
            [A,U,P,T,TA,TU,TP,AZT] = tube(M2,DT2,NBPTS,ELASTA,NB,SECTI,RAP,RO,A,U,P,T,C1,C2,C3,NBR,AZT);
            if nbif>0
                [A,U,P,T] = bif(DT2,PSING,nbif,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB);
            end
            if njonc>0
                [A,U,P,T] = jonc(DT2,PSING,nbif+1,njonc+nbif,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB);
            end
            if nelarg>0
                [A,U,P,T] = elarge(DT2,PSING,nbif+njonc+1,nbif+njonc+nelarg,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB);
            end
            if special==1
                [A,U,P,T] = big_jonction(DT2,BW,NBPTS,RAPI,ELASTA,SECTI,A,U,P,T,TA,TU,C1,C2,C3,NB);
            end
            [A,U,P,T,PE1,PIE,VVT,VVDT,QEI,QSI] = limit(N,DT2,RO,NBR,PRE,entree,sortie,fuite,A,U,P,T,ELASTA,SECTI,...
                TA,TU,TP,RAP,RAPI,NB,NBPTS,AZT,VVT,M2,C1,C2,C3,QEI,QSI,PE1,PIE,T_period,NPRE-1,PSORTIE);

            % Save into output files
            if (M==N)
                dis = sprintf('saving output... iteration N° %i',N);
                disp(dis);
                M=M+Msave;

                for br=1:NBR
                    i = 1;
                    DEB = A(br,i)*U(br,i);
                    fprintf(output_inlet,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
                end            
                for br=1:NBR
                    i = 4;
                    DEB = A(br,i)*U(br,i);
                    fprintf(output_pt4,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
                end
                for br=1:NBR
                    i = NBPTS(br)+1;
                    DEB = A(br,i)*U(br,i);
                    fprintf(output_outlet,'%12.5f  %12.5f  %12.5f  %12.5f  %12.5f\n',A(br,i),U(br,i),P(br,i),DEB,T(br,i));
                end    
            end
        end
        % -------------
        % END TIME LOOP
        % -------------

        IE=NBPTS(NBR)+1;

        QSN=sum(A(sortie',IE).*U(sortie',IE));
        if fuite>0
            QSN=QSN+A(fuite,1)*U(fuite,1);
        end
        QS=(QSI+0.5*(QS0-QSN))/N;
        QEN=sum(A(entree',1).*U(entree',1));
        QE=(QEI+0.5*(QE0-QEN))/N;
        dis = sprintf('QEI = %12.5f; QSI = %12.5f; VVT= %12.5f',QEI,QSI,VVT);
        disp(dis);

        VVM=(VVT-0.5*VVDT)/N;
        CQM=QS-QE+VVM;
        dis = sprintf('VVDT = %12.5f; VVM = %12.5f; QE= %12.5f; QS= %12.5f; CQM= %12.5f',VVDT,VVM,QE,QS,CQM);
        disp(dis);

        VOLT=0;
        for br=1:NBR
             IL=NBPTS(br)+1;
             VOLTI=0.5*LONG(br)*(A(br,1)+A(br,IL));
             VOLT=VOLT+VOLTI;
             dis = sprintf('branche N° %i; volume=%12.5f',br,VOLTI);
             disp(dis);
        end
        dis = sprintf('Volume total = %12.5f',VOLT);
        disp(dis);

    end

    fclose(output_inlet);
    fclose(output_pt4);
    fclose(output_outlet);

    % -----------------------------
    % Visualisations of the results
    % -----------------------------
    visua = input('Do you want to visualize the results? [y]es or [n]o   ','s');
    if ((visua=='y')|(visua=='Y'))
        VISUAL(FNAME);
    end
    
end


% -------------------------------------------------------------------------
%                       OTHER FUNCTIONS
% -------------------------------------------------------------------------

% -------------
% Vessels data
% -------------
function [NBR,DAD,SECTI,ELASTA,LONG,NB,NBPTS,entree,sortie,T_period,PS,NCycle,DT1,Msave] = data(datafile)

    fscanf(datafile,'%s',5);
    T_period = fscanf(datafile,'%f',1);
    PS = fscanf(datafile,'%f',1);
    NCycle = fscanf(datafile,'%i',1);
    DT1 = fscanf(datafile,'%f',1);
    Msave = fscanf(datafile,'%i',1);
    fscanf(datafile,'%s',4);
    NBR = fscanf(datafile,'%i');
    fscanf(datafile,'%s',3);
    howmany=fscanf(datafile,'%i',1);
    sortie=fscanf(datafile,'%i',howmany);
    fscanf(datafile,'%s',3);
    DAD = fscanf(datafile,'%i',NBR);
    fscanf(datafile,'%s',3);
    SECTI = fscanf(datafile,'%f',NBR);
    fscanf(datafile,'%s',3);
    ELASTA = fscanf(datafile,'%f',NBR);
    unitelasta = fscanf(datafile,'%f',1);
    ELASTA = unitelasta*ELASTA;
    fscanf(datafile,'%s',3);
    LONG = fscanf(datafile,'%f',NBR);
    fscanf(datafile,'%s',5);
    NB = fscanf(datafile,'%i',NBR);
    fscanf(datafile,'%s',4);
    NBPTS = fscanf(datafile,'%i',NBR);

    entree=find(DAD==0);
    
    DAD=DAD';
    SECTI=SECTI';
    ELASTA=ELASTA';
    LONG=LONG';
    NB=NB';
    NBPTS=NBPTS';

end

% ---------------------------
% INPUT MEASURE PRESSURE
% (to be adapted to measures)
% ---------------------------
function [PRE,NPRE] = entry()
    
    % Measurments points
    t=linspace(0,2*pi,86);
    PRE = 85.0+0.5*sin(t);
    NPRE = length(PRE);

end

% --------------
% Connectivities
% --------------
function [PSING,nbif,njonc,nelarg,fuite,special,BW] = connectivities(datafile)

    fscanf(datafile,'%s',7);
    nbif = fscanf(datafile,'%i',1);
    if nbif>0
        for i=1:nbif
            scan = fscanf(datafile,'%i',3);
            for j=1:3
               PSING(i,j) = scan(j);
           end
        end
    end
    
    fscanf(datafile,'%s',2);
    njonc = fscanf(datafile,'%i',1);
    if njonc>0
        for i=1:njonc
            scan = fscanf(datafile,'%i',3);
            for j=1:3
               PSING(nbif+i,j) = scan(j);
           end
        end
    end
    
    fscanf(datafile,'%s',2);
    nelarg = fscanf(datafile,'%i',1);
    if nelarg>0
        for i=1:nelarg
            scan = fscanf(datafile,'%i',3);
            for j=1:3
               PSING(nbif+njonc+i,j) = scan(j);
           end
        end
    end
    fscanf(datafile,'%s',3);
    special = fscanf(datafile,'%i',1);
    if special==1
        fscanf(datafile,'%s',2);
        howmany = fscanf(datafile,'%i',1);
        BW = fscanf(datafile,'%i',howmany);
    else
        BW = 0;
    end
    
    fuite = testfuite(PSING,nbif);
    
end


%    **********************************
%         NETWORK INITIALISATION
%       (to be adapted to the network)
%    **********************************
function [A,U,P,T] = init(entree,sortie,DAD,ELASTA,SECTI,LONG,NB,NBPTS,ETA,NBR,Z)
      
    QIN(1) = 3.0;
    QIN(2) = 3.0;
    QIN(3) = 0.0;
    QIN(4) = 3.0;    
    QIN(5) = 3.0;
    QIN(6) = 3.0;
    QIN(7) = 3.0;
    PIN = 85.0;
    POUT= 84.0;
    
    mmHg=1334;

    %   DONNEE PRESSION D'ENTREE
    PE=PIN*mmHg;

    %  CALCUL PRESSION SORTIE DES ENTREE
    PIN1(entree')=PE;
    AIN(entree')=(PIN1(entree')./ELASTA(entree')+1).*SECTI(entree');
    PIN2(entree')=PIN1(entree')-Z*QIN(entree').*LONG(entree').*NB(entree')./(AIN(entree').^2);
    
    % CALCUL PRESSION ENTREE /SORTIE DANS LE RESTE DU RESEAU    
    [sd,isort] = sort(DAD);    
    for nbr=entree(end)+1:NBR
        br=isort(nbr);
        PIN1(br)=PIN2(DAD(br));
        AIN=(PIN1(br)/ELASTA(br)+1)*SECTI(br);
        PIN2(br)=PIN1(br)-Z*QIN(br)*LONG(br)*NB(br)/(AIN^2);
    end

    % DONNEE PRESSION SORTIE RESEAU
    PSJ=POUT*mmHg;
    PIN2(sortie')=PSJ;
    PI1=PIN1/mmHg;
    PI2=PIN2/mmHg;

    % OUTPUT FILE: PI1,PI2,QIN
    output_file = fopen('PI12_QIN.dat','w');
    for i=1:NBR
        fprintf(output_file,'%12.5f  %12.5f  %12.5f\n',PI1(i),PI2(i),QIN(i));
    end
    fclose(output_file);


    % CALCUL PRESSION INITIAL EN TOUS LES PTS INTERNES DU RESEAU
    HG=1;
    for br =1:NBR
        P(br,1) = PIN1(br)*HG;
        P(br,NBPTS(br)+1)=PIN2(br)*HG;
    end
    for br=1:NBR
        IN = NBPTS(br)+1;
        DPI = (P(br,1)-P(br,IN))/(IN-1);
        for i=2:IN-1
            P(br,i)=P(br,i-1)-DPI;
        end
    end

    for br=1:NBR
        IN=NBPTS(br)+1;
        for i=1:IN
            A(br,i)=(P(br,i)/ELASTA(br)+1)*SECTI(br);
            U(br,i)=QIN(br)/A(br,i);
            T(br,i)=-4*ETA*sqrt(pi*NB(br)/A(br,i))*U(br,i);
        end
    end

end

% -----------------------------------
% test to find if there's a leak tube
% -----------------------------------
function [fuite] = testfuite(PSING,nbif)

    fuite=0;
    for i=1:nbif
        if PSING(i,2)==PSING(i,3)
            fuite=PSING(i,2);
        end
    end
    
end


%    *******************************
%          CALCUL PTS INTERNES
%    *******************************
function [A,U,P,T,TA,TU,TP,AZT] = tube(M1,DT1,NBPTS,ELASTA,NB,SECTI,RAP,RO,A,U,P,T,C1,C2,C3,NBR,AZT)
      
    for br=1:NBR

        k=br;
        M=M1;
        DT=DT1;
        I3=NBPTS(k);
        I1=1;
        I2=2;
        I4=I3+1;
        EL=ELASTA(k);
        NBB=NB(k);
        RNB=sqrt(NBB);
        A0=SECTI(k);
        R=RAP(k);
        R2=R/2;
        M3=M/2;
        Z1=DT*sqrt(pi)*RNB/RO;


    %   CALCUL 1ER DEMI-PAS

        for I=I1:I3

          AA=A(k,I);
          UU=U(k,I);
          PP=P(k,I);
          AA1=A(k,I+1);
          UU1=U(k,I+1);
          PP1=P(k,I+1);
          TT=T(k,I);
          TT1=T(k,I+1);

          AM=0.5*(AA+AA1);
          UM=0.5*(UU+UU1);
          TM=0.5*(TT+TT1);
          RAM=sqrt(AM);


          A1Z=R2*(AA1*UU1-AA*UU);
          A1(I+1)=(AA1+AA)/2.-A1Z;
          Q1=UU1*UU1/2+PP1/RO;
          Q=UU*UU/2+PP/RO;
          U1(I+1)=UM-R2*(Q1-Q)+Z1*TM/RAM;

          B2=C2*RNB*UM/RAM;
          B3=(C3*RAM/RNB)*((U1(I+1)-UM)/DT/2);

          T1(I+1)=TM-(0.5*(DT*NBB)/(C1*AM))*(TM+B2+B3);
          P1(I+1)=EL*(A1(I+1)/A0-1);
        end

    %   CALCUL 2EME DEMI-PAS
        I=I1;
        TA(k,1)=A(k,1);
        TU(k,1)=U(k,1);
        TP(k,1)=P(k,1);
        TA(k,2)=A(k,2);
        TU(k,2)=U(k,2);
        TP(k,2)=P(k,2);
        I=I3;
        TA(k,3)=A(k,I3);
        TU(k,3)=U(k,I3);
        TP(k,3)=P(k,I3);
        I=I4;
        TA(k,4)=A(k,I4);
        TU(k,4)=U(k,I4);
        TP(k,4)=P(k,I4);

        for I=I2:I3

          AA=A1(I);
          UU=U1(I);
          PP=P1(I);
          AA1=A1(I+1);
          UU1=U1(I+1);
          PP1=P1(I+1);
          TT=T1(I);
          TT1=T1(I+1);

          AM=0.5*(AA+AA1);
          UM=0.5*(UU+UU1);
          TM=0.5*(TT+TT1);
          RAM=sqrt(AM);

          AZ=-R*(AA1*UU1-AA*UU);

          AZT=AZT+AZ/R;

          A(k,I)=A(k,I)+AZ;
          S1=UU1*UU1/2+PP1/RO;
          S=UU*UU/2+PP/RO;

          UN=U(k,I);
          U(k,I)=UN-R*(S1-S)+2*Z1*TM/RAM;

          B2=C2*RNB*UM/RAM;
          B3=(C3*RAM/RNB)*((U(k,I)-UN)/DT);

          T(k,I)=T(k,I)-((DT*NBB)/(C1*AM))*(TM+B2+B3);
          P(k,I)=EL*(A(k,I)/A0-1);

        end

    end

end


%    ***************************************
%         TRAITEMENT DES BIFURCATIONS
%    ****************************************
function [A,U,P,T] = bif(DT2,PSING,nbif,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB)
      
    DT=DT2;
    for PS=1:nbif
        BP=PSING(PS,1);
        BS=PSING(PS,2);
        BF=PSING(PS,3);
        I0=NBPTS(BP);
        I4=I0+1;
        I1=1;
        I2=2;
        I1F=1;
        I2F=2;
        R1=RAPI(BP)/2;
        R2=RAPI(BS)/2;
        R3=RAPI(BF)/2;
        A0P=SECTI(BP);
        A0S=SECTI(BS);
        A0F=SECTI(BF);
        EL1=ELASTA(BP);
        EL2=ELASTA(BS);
        EL3=ELASTA(BF);
        KCS=(EL1*A0S)/(EL2*A0P);
        KCF=(EL1*A0F)/(EL3*A0P);

        T1=A(BP,I0)*U(BP,I0);
        T2=A(BS,I2)*U(BS,I2)+A(BF,I2F)*U(BF,I2F);
        T3=R1*(A(BP,I0)-TA(BP,3));
        T4=R2*(A(BS,I2)-TA(BS,2))+R3*(A(BF,I2F)-TA(BF,I2F));
        T5=R1+R2*KCS+R3*KCF;
        T6=R1-R2*KCS-R3*KCF;

        AC=TA(BP,4)+(T1-T2-T3-T4)/T5;
        UC=(T1+T2+T4-T3-T6*(AC-TA(BP,4)))/(2*AC);
        PC=EL1*(AC/A0P-1);

        if (BS==BF)
            ACS=(1+PC/EL2)*A0S;
            UCS=AC*UC/(2*ACS);
            UCF=UCS;
            ACF=ACS;
        else
            ACS=(1+PC/EL2)*A0S;
            UCS=(R2*(A(BS,I2)-TA(BS,2)+ACS-TA(BS,1))+A(BS,I2)*U(BS,I2))/ACS;
            ACF=(1+PC/EL3)*A0F;
            UCF=(AC*UC-ACS*UCS)/ACF;
        end

        [T] = stress(4,AC,UC,DT,BP,NB,NBPTS,C1,C2,C3,T,TU);
        [T] = stress(1,ACS,UCS,DT,BS,NB,NBPTS,C1,C2,C3,T,TU);
        [T] = stress(1,ACF,UCF,DT,BF,NB,NBPTS,C1,C2,C3,T,TU);

        A(BP,I4)=AC;
        U(BP,I4)=UC;
        P(BP,I4)=PC;
        A(BS,I1)=ACS;
        U(BS,I1)=UCS;
        P(BS,I1)=PC;
        P(BF,I1F)=PC;
        A(BF,I1F)=ACF;
        U(BF,I1F)=UCF;
    end
    
end


%**************************************************************
%     CALCUL WALL SHEAR STRESS ELARGE,BIF,JONC,BIG_JONCTION
%**************************************************************
function [T] = stress(J1,AN1,UN1,DT,BR,NB,NBPTS,C1,C2,C3,T,TU)

    NBB=NB(BR);

    if (J1==1)
        I=J1;
    elseif (J1>1)
        I=NBPTS(BR)+1;
    else
        return;
    end
    B0=(C1*AN1)/(NBB*DT);
    B1=B0+1;
    B2=C2*sqrt(NBB/AN1)*UN1;
    B3=C3*sqrt(AN1/NBB)*(UN1-TU(BR,J1))/DT;
    T(BR,I)=(B0*T(BR,I)-B2-B3)/B1;
    %T(BR,I)=-C2*SQRT(NBB/TA(BR,J1))*TU(BR,J1)

end

%    ***************************************
%          TRAITEMENT DES JONCTIONS
%    ****************************************
function [A,U,P,T] = jonc(DT2,PSING,debut,fin,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB)
      
    DT=DT2;

    for PS=debut:fin

        BP=PSING(PS,1);
        BS=PSING(PS,2);
        BF=PSING(PS,3);
        I1=1;
        I2=2;
        I0F=NBPTS(BF);
        IF=I0F+1;
        I0S=NBPTS(BS);
        IS=I0S+1;
        R1=RAPI(BP)/2;
        R2=RAPI(BS)/2;
        R3=RAPI(BF)/2;
        A0P=SECTI(BP);
        A0S=SECTI(BS);
        A0F=SECTI(BF);
        EL1=ELASTA(BP);
        EL2=ELASTA(BS);
        EL3=ELASTA(BF);
        KCS=(EL1*A0S)/(EL2*A0P);
        KCF=(EL1*A0F)/(EL3*A0P);

        T1=A(BP,I2)*U(BP,I2);
        T2=A(BS,I0S)*U(BS,I0S)+A(BF,I0F)*U(BF,I0F);
        T3=R1*(A(BP,I2)-TA(BP,2));
        T4=R2*(A(BS,I0S)-TA(BS,3))+R3*(A(BF,I0F)-TA(BF,3));
        T5=R1+R2*KCS+R3*KCF;
        T6=R1-R2*KCS-R3*KCF;
        AC=(T2-T1-T3-T4)/T5+TA(BP,1);
        UC=(T1+T2+T3-T4+T6*(AC-TA(BP,1)))/(2*AC);
        PC=EL1*(AC/A0P-1);
        ACS=(1+PC/EL2)*A0S;

        if (SECTI(BS)==SECTI(BF))
            UCS=AC*UC/(2*ACS);
            UCF=UCS;
            ACF=ACS;
        else
            UCS=(A(BS,I0S)*U(BS,I0S)-R2*(A(BS,I0S)-TA(BS,3)+ACS-TA(BS,4)))/ACS;
            ACF=(1+PC/EL3)*A0F;
            UCF=(AC*UC-ACS*UCS)/ACF;
        end

        [T] = stress(1,AC,UC,DT,BP,NB,NBPTS,C1,C2,C3,T,TU);
        [T] = stress(4,ACS,UCS,DT,BS,NB,NBPTS,C1,C2,C3,T,TU);
        [T] = stress(4,ACF,UCF,DT,BF,NB,NBPTS,C1,C2,C3,T,TU);

        A(BP,I1)=AC;
        U(BP,I1)=UC;
        P(BP,I1)=PC;
        A(BS,IS)=ACS;
        U(BS,IS)=UCS;
        P(BS,IS)=PC;
        A(BF,IF)=ACF;
        U(BF,IF)=UCF;
        P(BF,IF)=PC;

    end
end

%    ********************************
%     TRAITEMENT DES ELARGISSEMENTS
%        ET DES RETRECISSEMENTS
%    ********************************
function [A,U,P,T] = elarge(DT2,PSING,debut,fin,NBPTS,RAPI,SECTI,ELASTA,A,U,P,T,TA,TU,C1,C2,C3,NB)
      
    for PS=debut:fin

        DT=DT2;
        BP=PSING(PS,1);
        BS=PSING(PS,2);
        I0=NBPTS(BP);
        I1=1;
        I2=2;
        I4=I0+1;
        R1=RAPI(BP)/2;
        R2=RAPI(BS)/2;

        A0P=SECTI(BP);
        A0S=SECTI(BS);
        EL1=ELASTA(BP);
        EL2=ELASTA(BS);
        KCS=(EL1*A0S)/(EL2*A0P);
        T1=A(BP,I0)*U(BP,I0);
        T2=A(BS,I2)*U(BS,I2);
        T3=R1*(A(BP,I0)-TA(BP,3));
        T4=R2*(A(BS,2)-TA(BS,2));
        T5=R1+R2*KCS;
        T6=R1-R2*KCS;
        AC=TA(BP,4)+(T1-T2-T3-T4)/T5;
        UC=(T1+T2+T4-T3-T6*(AC-TA(BP,4)))/(2*AC);
        PC=EL1*(AC/A0P-1);

        [T] = stress(4,AC,UC,DT,BP,NB,NBPTS,C1,C2,C3,T,TU);

        ACP=(1+PC/EL2)*A0S;
        UCP=UC*AC/ACP;

        [T] = stress(1,ACP,UCP,DT,BS,NB,NBPTS,C1,C2,C3,T,TU);

        A(BP,I4)=AC;
        U(BP,I4)=UC;
        P(BP,I4)=PC;
        A(BS,I1)=ACP;
        U(BS,I1)=UCP;
        P(BS,I1)=PC;

    end

end

%    ***************************
%     TRAITEMENT D'UNE JONCTION 
%     DE PLUSIEURS VAISSEAUX
%    ****************************
function [A,U,P,T] = big_jonction(DT2,BW,NBPTS,RAPI,ELASTA,SECTI,A,U,P,T,TA,TU,C1,C2,C3,NB)

    DT=DT2;
    NPTW=length(BW);

    for i=1:NPTW
        B(i)=BW(i);
        I3(i)=NBPTS(B(i));
        I4(i)=I3(i)+1;
        R(i)=RAPI(B(i))/2;
        EL(i)=ELASTA(B(i));
        A0(i)=SECTI(B(i));
    end
    I1(NPTW)=1;
    I2(NPTW)=I1(NPTW)+1;
    
    for i=1:NPTW-1
        KC(i)=(EL(NPTW)*A0(i))/(EL(i)*A0(NPTW));
    end

    T1=A(B(NPTW),I2(NPTW))*U(B(NPTW),I2(NPTW));
    T2=0.0;
    for i=1:NPTW-1
        T2=T2+A(B(i),I3(i))*U(B(i),I3(i));
    end
    T3=R(NPTW)*(A(B(NPTW),I2(NPTW))-TA(B(NPTW),2));
    T4=0.0;
    for i=1:NPTW-1
        T4=T4+R(i)*(A(B(i),I3(i))-TA(B(i),3));
    end
    T5=R(NPTW)+sum(R(1:NPTW-1).*KC(1:NPTW-1));
    T6=R(NPTW)-sum(R(1:NPTW-1).*KC(1:NPTW-1));

    AC(NPTW)=(T2-T1-T3-T4)/T5+TA(B(NPTW),1);
    UC(NPTW)=(T1+T2+T3-T4+T6*(AC-TA(B(NPTW),1)))/(2*AC);
    PC=EL(NPTW)*(AC(NPTW)/A0(NPTW)-1);
    
    for i=1:NPTW-1
        AC(i)=(1+PC/EL(i))*A0(i);
        UC(i)=(A(B(i),I3(i))*U(B(i),I3(i))-R(i)*(A(B(i),I3(i))-TA(B(i),3)+AC(i)-TA(B(i),4)))/AC(i);
        [T] = stress(4,AC(i),UC(i),DT,B(i),NB,NBPTS,C1,C2,C3,T,TU);
    end
    [T] = stress(1,AC(NPTW),UC(NPTW),DT,B(NPTW),NB,NBPTS,C1,C2,C3,T,TU);

    for i=1:NPTW-1
        A(B(i),I4(i))=AC(i);
        U(B(i),I4(i))=UC(i);
        P(B(i),I4(i))=PC;
    end
    A(B(NPTW),I1(NPTW))=AC(NPTW);
    U(B(NPTW),I1(NPTW))=UC(NPTW);
    P(B(NPTW),I1(NPTW))=PC;

end


%    ****************************************
%     CONDITIONS AUX LIMITES ENTREES_SORTIES
%    *****************************************
function [A,U,P,T,PE1,PIE,VVT,VVDT,QEI,QSI] = limit(N,DT2,RO,NBR,PRE,entree,sortie,fuite,A,U,P,T,ELASTA,SECTI,...
    TA,TU,TP,RAP,RAPI,NB,NBPTS,AZT,VVT,M2,C1,C2,C3,QEI,QSI,PE1,PIE,T_period,DATAT,PSORTIE)

    DT=DT2;
    Z1=DT*sqrt(pi)/RO;
    N1=N;
    e=ceil(T_period/DT2)-fix(T_period/DT2);
    DELTAT=fix(T_period/DT2+e)/DATAT;
    I=1+fix((N1-1)/DELTAT);
    X=N1-DELTAT*(I-1);

    if (X<=1)
        PE1=PRE(I);
        PE2=PRE(I+1);
        PIE=(PE2-PE1)/DELTAT;
    end

    I1=NBPTS(NBR)+1;

    PE=(PE1+X*PIE)*1334;

    for k=1:length(entree)
        BR=entree(k);
        P(BR,1)=PE;
        A(BR,1)=(PE/ELASTA(BR)+1)*SECTI(BR);

        S1=TU(BR,1)*TU(BR,1)/2.+TP(BR,1)/RO;
        S2=TU(BR,2)*TU(BR,2)/2.+TP(BR,2)/RO;

        U(BR,1)=TU(BR,1)+RAP(BR)*(S1-S2)-M2*TU(BR,1)/TA(BR,1);
        U(BR,1)=TU(BR,1)-RAP(BR)*(S2-S1)+2*Z1*T(BR,1)/sqrt(TA(BR,1));

        UN1=U(BR,1);
        [T] = stressl(1,UN1,DT,BR,NB,NBPTS,T,TA,TU,C1,C2,C3);
    end

    for k=1:length(sortie)
        BR=sortie(k);
        IL=NBPTS(BR)+1;
        S1=TU(BR,3)*TU(BR,3)/2.+TP(BR,3)/RO;
        S2=TU(BR,4)*TU(BR,4)/2.+TP(BR,4)/RO;

        U(BR,IL)=TU(BR,4)-RAP(BR)*(S2-S1)+2*Z1*T(BR,IL)/sqrt(TA(BR,4));

        UN1=U(BR,IL);
        [T] = stressl(4,UN1,DT,BR,NB,NBPTS,T,TA,TU,C1,C2,C3);
        P(BR,IL)=PSORTIE*1334;
        A(BR,IL)=SECTI(BR)*(P(BR,IL)/ELASTA(BR)+1);
    end

    VAL=0;
    for BR=1:NBR
          IL=NBPTS(BR)+1;
          VAL=VAL+0.5*RAPI(BR)*(A(BR,1)-TA(BR,1)+A(BR,IL)-TA(BR,4));
    end
    VVDT=VAL+AZT;
    VVT=VVT+VVDT;

    QLE=sum(A(entree',1).*U(entree',1));

    IL=NBPTS(NBR)+1;
    QEI=QEI+QLE;
    QSI=QSI+sum(A(sortie',IL).*U(sortie',IL));
    if fuite>0
        QSI=QSI+A(fuite,1)*U(fuite,1);
    end
end


%    ****************************************
%     CALCUL WALL SHEAR STRESS IN/OUT DU RESEAU
%    *****************************************
function [T] = stressl(J1,UN1,DT,BR,NB,NBPTS,T,TA,TU,C1,C2,C3)

    NBB=NB(BR);

    if (J1==1)
        I=J1;
    elseif (J1>1)
        I=NBPTS(BR)+1;
    else
        return;
    end

    B1=(DT*NBB)/(C1*TA(BR,J1));
    B2=C2*sqrt(NBB/TA(BR,J1))*TU(BR,J1);
    B3=C3*sqrt(TA(BR,J1)/NBB)*(UN1-TU(BR,J1))/DT;
    T(BR,I)=T(BR,I)-B1*(T(BR,I)+B2+B3);

end