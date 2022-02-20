%**************************************************************************
%
%           Solution analytique de la propagation
%            des ondes en conduites déformables
%
%       (from Chirazi, Khalid Naciri et Ouazzani Touhami, 1998)
% A(x=0,t) = A0+a0*H(t); A(x=L,t) = A0+a1*B(t)
% A(x,t=0) = Ai(x); Q(x,t=0) = Qi(x)
%
% (P. Cathalifaud, 2022)
%**************************************************************************

function analytical_solution()

    global L A0 T0 V0 a0 a1;
    
    % Longueur du tube
    L = 25;
    % Section au repos
    A0 = 0.14;
 
    % Perturbation de section en entrée et sortie
    a0 = 0.014;
    a1 = 0.014;
    
    % Temps caractéristique du signal de perturbation
    T0 = 0.08;
    
    % Coefficient de souplesse du tube
    XI = 55*1.0E-7; % EL = XI/A0
    
    % Densité du fluide
    RHO = 1.05;
    
    % Vitesse de propagation et temps de parcours du tube des ondes
    V0 = sqrt(A0/(RHO*XI));
    temps_ondes = L/V0;
    
    fprintf('%f m/s: Vitesse des ondes \n',V0);
    fprintf('%f s: Temps de parcours du tube par les ondes \n',temps_ondes);
    
    t = linspace(0,0.8,1000);
    
    % Position de visualisation dans le tube
    z = L/4;
    
    % Solutions analytiques
    for i=1:length(t)
        A(i) = reflected(z+V0*t(i)) + incidente(z-V0*t(i));
        f(i) = reflected(z+V0*t(i));
        g(i) = incidente(z-V0*t(i));
    end
    
    % Coefficient pour changer d'unité (cm2 --> mm2)
    mm2 = 100;
    
    % Figures
    figure;plot(t,A*mm2,'k','linewidth',2);
    grid on;title('Section totale (mm^2)','FontSize',30);xlabel('time (s)','FontSize',20);
    figure;plot(t,f*mm2,'k','linewidth',2);
    grid on;title('Section réfléchie (mm^2)','FontSize',30);xlabel('time (s)','FontSize',20);
    figure;plot(t,g*mm2,'k','linewidth',2);
    grid on;title('Section incidente (mm^2)','FontSize',30);xlabel('time (s)','FontSize',20);
        
end

function sortie = H(t)

    global T0;
        
    if t<=T0
        sortie = sin(2*pi*t/T0);
    else
        sortie = 0.0;
    end
    
end

function sortie = B(t)

    sortie = 0.0;
    
end

function sortie = Qi(x)
    
    sortie = 0.0;
    
end

function sortie = Ai(x)

    global A0;

    sortie = A0;
    
end

function f = reflected(s)

    global L V0 A0 a0 a1;
    
    n = 0;
    while s>n*L
        n = n+1;
    end
    
    if mod(n,2)~=0
        k = (n-1)/2;
        terme1=0;
        for m=0:k-1
            terme1 = terme1+a1*B((s-(2*m+1)*L)/V0);
        end
        terme2=0;
        for m=1:k
            terme2 = terme2+a0*H((s-2*m*L)/V0);
        end
        f = -1/(2*V0)*Qi(s-2*k*L)+Ai(s-2*k*L)/2+terme1-terme2-A0/2;
    else
        k = (n-2)/2;
        terme1=0;
        for m=0:k
            terme1 = terme1+a1*B((s-(2*m+1)*L)/V0);
        end
        terme2=0;
        for m=1:k
            terme2 = terme2+a0*H((s-2*m*L)/V0);
        end
        f = A0-1/(2*V0)*Qi((2*k+2)*L-s)-Ai((2*k+2)*L-s)/2+terme1-terme2-A0/2;
    end
    
end

function g = incidente(s)

    global A0 a0 V0;

    if s<=0
        g = A0-reflected(-s)+a0*H(-s/V0);
    else
        g = Ai(s)-reflected(s);
    end
    
end
