%************************************************
%           Bifurcation 0D avec ballons
%
%                          O <-- Ballon V
%                         /
%                        /
%   ENTREE -->  ________/
%                       \
%                        \
%                         \
%                          O <-- Ballon E
%
% P. Cathalifaud
%************************************************

%% MAIN FUNCTION
function BIFURCATION0D()

     cv = 2;
     ce = 1;
     Le = 100;
     Lv = Le;
     rho = 1.0E-03;
     mu = 1.0E-03;
     Re = 5;
     Rv = 1;
     Ae = pi*Re^2;
     Av = pi*Rv^2;
     Qmax = 10;
     T = 1;
     omega = 2*pi/T;


    cas = input('Type approx pour le terme visqueux ? steady/Poiseuille (1) or unsteady/Womersley (2) ' );
    if (cas==1)
        % Poiseuille
        sigmae=8*pi*mu*Le/(Ae^2);
        sigmav=8*pi*mu*Lv/(Av^2);
    else
        % Approximation des sigmas (termes visqueux)
        % Unsteady
        worme = Re*sqrt(omega);
        wormv = Rv*sqrt(omega);
        F10e = 2*besselj(1,i^(3/2)*worme)/(i^(3/2)*worme*besselj(0,i^(3/2)*worme));
        F10v = 2*besselj(1,i^(3/2)*wormv)/(i^(3/2)*wormv*besselj(0,i^(3/2)*wormv));

        sigmae=i*rho*omega*(F10e./(1-F10e))*Le/Ae;
        sigmav=i*rho*omega*(F10v./(1-F10v))*Lv/Av;
    end


    alpha1=rho*(Le/Ae+Lv/Av);
    alpha2=sigmae+sigmav;
    alpha3=1/ce+1/cv;

    % Termes de forçage
    Kv1=Qmax*(1/ce-rho*Le*omega^2/Ae);
    Kv2=Qmax*sigmae*omega;
    Ke1=Qmax*(1/cv-rho*Lv*omega^2/Av);
    Ke2=Qmax*sigmav*omega;

    dt = 0.001;
    t=[18:dt:20];

    if (cas==1)
        % Cas forçage réel: Qmax*sin(omega*t)
        M = [alpha3-alpha1*omega^2,-alpha2*omega;alpha2*omega,alpha3-alpha1*omega^2];
        [AvBv] = inv(M)*[Kv1;Kv2];
        [AeBe] = inv(M)*[Ke1;Ke2];
        Qe = AeBe(1)*sin(omega*t)+AeBe(2)*cos(omega*t);
        Qv = AvBv(1)*sin(omega*t)+AvBv(2)*cos(omega*t);
    else
        % Cas forçage imaginaire: Qmax*exp(i*omega*t)
        XIe = (Ke1+i*Ke2)/(alpha3-omega^2*alpha1+i*omega*alpha2);
        Qe_tilde = sqrt( real(XIe)^2+imag(XIe)^2 );
        Phie = atan2( imag(XIe),real(XIe) );
        Qe = Qe_tilde*sin(omega*t+Phie);
        XIv = (Kv1+i*Kv2)/(alpha3-omega^2*alpha1+i*omega*alpha2);
        Qv_tilde = sqrt( real(XIv)^2+imag(XIv)^2 );
        Phiv = atan2( imag(XIv),real(XIv) );
        Qv = Qv_tilde*sin(omega*t+Phiv);
    end

    Qsomme = Qe+Qv;
    Q = Qmax*sin(omega*t);
    figure;plot(t,Qe,'r',t,Qv,'g',t,Q,'k',t,Qsomme,'b--');grid on;

    % Calcul du stroke volume
    Se = 0;
    Sv = 0;
    S = 0;
    for ind=1:length(Qe)-1
        Se = Se + (abs(Qe(ind)+Qe(ind+1)))/2*dt;
        Sv = Sv + (abs(Qv(ind)+Qv(ind+1)))/2*dt;
        S = S + (abs(Q(ind)+Q(ind+1)))/2*dt;
    end

    fprintf('STROKE VOLUMES \n');
    fprintf('SAS VENTRICULE');
    [Se/(Se+Sv) Sv/(Se+Sv)]

    clear omega;
    clear t;

    omega = [0:0.1:50];


    for ind1=1:length(omega)
        clear t;
        T = 2*pi/omega(ind1);
        dt = T/100;
        t=[0:dt:T];

        worme = Re*sqrt(omega(ind1));
        wormv = Rv*sqrt(omega(ind1));
        F10e = 2*besselj(1,i^(3/2)*worme)/(i^(3/2)*worme*besselj(0,i^(3/2)*worme));
        F10v = 2*besselj(1,i^(3/2)*wormv)/(i^(3/2)*wormv*besselj(0,i^(3/2)*wormv));

        sigmae=i*rho*omega(ind1)*(F10e./(1-F10e))*Le/Ae;
        sigmav=i*rho*omega(ind1)*(F10v./(1-F10v))*Lv/Av;

        alpha1=rho*(Le/Ae+Lv/Av);
        alpha2=sigmae+sigmav;
        alpha3=1/ce+1/cv;

        % Termes de forçage
        Kv1=Qmax*(1/ce-rho*Le*omega(ind1)^2/Ae);
        Kv2=Qmax*sigmae*omega(ind1);
        Ke1=Qmax*(1/cv-rho*Lv*omega(ind1)^2/Av);
        Ke2=Qmax*sigmav*omega(ind1);

        XIe = (Ke1+i*Ke2)/(alpha3-omega(ind1)^2*alpha1+i*omega(ind1)*alpha2);
        Qe_tilde = sqrt( real(XIe)^2+imag(XIe)^2 );
        Phie = atan2( imag(XIe),real(XIe) );

        XIv = (Kv1+i*Kv2)/(alpha3-omega(ind1)^2*alpha1+i*omega(ind1)*alpha2);
        Qv_tilde = sqrt( real(XIv)^2+imag(XIv)^2 );
        Phiv = atan2( imag(XIv),real(XIv) );

        Qe = Qe_tilde*sin(omega(ind1)*t+Phie);
        Qv = Qv_tilde*sin(omega(ind1)*t+Phiv);


        Se = 0;
        Sv = 0;
        S = 0;
        for ind=1:length(Qe)-1
            Se = Se + (abs(Qe(ind)+Qe(ind+1)))/2*dt;
            Sv = Sv + (abs(Qv(ind)+Qv(ind+1)))/2*dt;
            S = S + (abs(Q(ind)+Q(ind+1)))/2*dt;
        end
        Ste(ind1) = Se;
        Stv(ind1) = Sv;
        St(ind1) = S;
    end

    Ste_adim = Ste./(Ste+Stv);
    Stv_adim = Stv./(Ste+Stv);

    figure;
    plot(omega/(2*pi)*60,Ste_adim,'r',omega/(2*pi)*60,Stv_adim,'g');
    grid on;
    xlabel('rythme cardiaque [bpm]');
    ylabel('Stroke volume [%]');

end
