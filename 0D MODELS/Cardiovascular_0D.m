function [VOI, STATES, ALGEBRAIC, CONSTANTS] = Cardiovascular_0D()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =34;
end
% There are a total of 14 entries in each of the rate and state variable arrays.
% There are a total of 88 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 10];

    % Set numerical accuracy options for ODE solver
    options = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 1);

    % Solve model with ODE solver
    [VOI, STATES] = ode15s(@(VOI, STATES)computeRates(VOI, STATES, CONSTANTS), tspan, INIT_STATES, options);

    % Compute algebraic variables
    [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS);
    ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI);

    % Plot state variables against variable of integration
    [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends();
    figure();
    plot(VOI, STATES);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES);
    set(l,'Interpreter','none');
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('t in component environment (second)');
    LEGEND_ALGEBRAIC(:,19) = strpad('Pi in component TempCDa (UnitP)');
    LEGEND_STATES(:,1) = strpad('Pi in component TempRLC (UnitP)');
    LEGEND_ALGEBRAIC(:,31) = strpad('Qo in component TempRC (UnitQ)');
    LEGEND_ALGEBRAIC(:,11) = strpad('Qo in component TempCDv (UnitQ)');
    LEGEND_ALGEBRAIC(:,20) = strpad('Pi in component TempCDa (UnitP)');
    LEGEND_STATES(:,2) = strpad('Pi in component TempRLC (UnitP)');
    LEGEND_ALGEBRAIC(:,32) = strpad('Qo in component TempRC (UnitQ)');
    LEGEND_ALGEBRAIC(:,12) = strpad('Qo in component TempCDv (UnitQ)');
    LEGEND_ALGEBRAIC(:,7) = strpad('Pi in component TempCDv (UnitP)');
    LEGEND_ALGEBRAIC(:,23) = strpad('Qo in component TempCDa (UnitQ)');
    LEGEND_CONSTANTS(:,1) = strpad('CVao in component ParaHeart (UnitCV)');
    LEGEND_ALGEBRAIC(:,5) = strpad('E in component EVentricle (UnitE)');
    LEGEND_STATES(:,3) = strpad('V in component TempCDv (UnitV)');
    LEGEND_CONSTANTS(:,2) = strpad('PlvIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,3) = strpad('VlvIni in component ParaHeart (UnitV)');
    LEGEND_ALGEBRAIC(:,9) = strpad('Tao in component TempCDv (dimensionless)');
    LEGEND_CONSTANTS(:,4) = strpad('Vlv0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,5) = strpad('CVmi in component ParaHeart (UnitCV)');
    LEGEND_ALGEBRAIC(:,17) = strpad('E in component EAtrium (UnitE)');
    LEGEND_STATES(:,4) = strpad('V in component TempCDa (UnitV)');
    LEGEND_CONSTANTS(:,6) = strpad('PlaIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,7) = strpad('VlaIni in component ParaHeart (UnitV)');
    LEGEND_ALGEBRAIC(:,21) = strpad('Tao in component TempCDa (dimensionless)');
    LEGEND_CONSTANTS(:,8) = strpad('Vla0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,9) = strpad('ElvMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,10) = strpad('ElvMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,11) = strpad('T in component ParaHeart (second)');
    LEGEND_CONSTANTS(:,12) = strpad('Ts1 in component ParaHeart (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('Ts2 in component ParaHeart (dimensionless)');
    LEGEND_ALGEBRAIC(:,1) = strpad('mt in component EVentricle (second)');
    LEGEND_ALGEBRAIC(:,3) = strpad('et in component EVentricle (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('ElaMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,15) = strpad('ElaMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,16) = strpad('Tpwb in component ParaHeart (dimensionless)');
    LEGEND_CONSTANTS(:,17) = strpad('Tpww in component ParaHeart (dimensionless)');
    LEGEND_ALGEBRAIC(:,13) = strpad('mt in component EAtrium (second)');
    LEGEND_ALGEBRAIC(:,15) = strpad('et in component EAtrium (dimensionless)');
    LEGEND_CONSTANTS(:,18) = strpad('EraMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,19) = strpad('EraMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,20) = strpad('PraIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,21) = strpad('VraIni in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,22) = strpad('ErvMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,23) = strpad('ErvMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,24) = strpad('PrvIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,25) = strpad('VrvIni in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,26) = strpad('CVpa in component ParaHeart (UnitCV)');
    LEGEND_CONSTANTS(:,27) = strpad('CVti in component ParaHeart (UnitCV)');
    LEGEND_CONSTANTS(:,28) = strpad('Vra0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,29) = strpad('Vrv0 in component ParaHeart (UnitV)');
    LEGEND_ALGEBRAIC(:,8) = strpad('Pi in component TempCDv (UnitP)');
    LEGEND_ALGEBRAIC(:,24) = strpad('Qo in component TempCDa (UnitQ)');
    LEGEND_CONSTANTS(:,30) = strpad('CVpa in component ParaHeart (UnitCV)');
    LEGEND_ALGEBRAIC(:,6) = strpad('E in component EVentricle (UnitE)');
    LEGEND_STATES(:,5) = strpad('V in component TempCDv (UnitV)');
    LEGEND_CONSTANTS(:,31) = strpad('PrvIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,32) = strpad('VrvIni in component ParaHeart (UnitV)');
    LEGEND_ALGEBRAIC(:,10) = strpad('Tao in component TempCDv (dimensionless)');
    LEGEND_CONSTANTS(:,33) = strpad('Vrv0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,34) = strpad('CVti in component ParaHeart (UnitCV)');
    LEGEND_ALGEBRAIC(:,18) = strpad('E in component EAtrium (UnitE)');
    LEGEND_STATES(:,6) = strpad('V in component TempCDa (UnitV)');
    LEGEND_CONSTANTS(:,35) = strpad('PraIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,36) = strpad('VraIni in component ParaHeart (UnitV)');
    LEGEND_ALGEBRAIC(:,22) = strpad('Tao in component TempCDa (dimensionless)');
    LEGEND_CONSTANTS(:,37) = strpad('Vra0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,38) = strpad('ErvMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,39) = strpad('ErvMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,40) = strpad('T in component ParaHeart (second)');
    LEGEND_CONSTANTS(:,41) = strpad('Ts1 in component ParaHeart (dimensionless)');
    LEGEND_CONSTANTS(:,42) = strpad('Ts2 in component ParaHeart (dimensionless)');
    LEGEND_ALGEBRAIC(:,2) = strpad('mt in component EVentricle (second)');
    LEGEND_ALGEBRAIC(:,4) = strpad('et in component EVentricle (dimensionless)');
    LEGEND_CONSTANTS(:,43) = strpad('EraMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,44) = strpad('EraMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,45) = strpad('Tpwb in component ParaHeart (dimensionless)');
    LEGEND_CONSTANTS(:,46) = strpad('Tpww in component ParaHeart (dimensionless)');
    LEGEND_ALGEBRAIC(:,14) = strpad('mt in component EAtrium (second)');
    LEGEND_ALGEBRAIC(:,16) = strpad('et in component EAtrium (dimensionless)');
    LEGEND_CONSTANTS(:,47) = strpad('ElaMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,48) = strpad('ElaMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,49) = strpad('PlaIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,50) = strpad('VlaIni in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,51) = strpad('ElvMax in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,52) = strpad('ElvMin in component ParaHeart (UnitE)');
    LEGEND_CONSTANTS(:,53) = strpad('PlvIni in component ParaHeart (UnitP)');
    LEGEND_CONSTANTS(:,54) = strpad('VlvIni in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,55) = strpad('CVao in component ParaHeart (UnitCV)');
    LEGEND_CONSTANTS(:,56) = strpad('CVmi in component ParaHeart (UnitCV)');
    LEGEND_CONSTANTS(:,57) = strpad('Vlv0 in component ParaHeart (UnitV)');
    LEGEND_CONSTANTS(:,58) = strpad('Vla0 in component ParaHeart (UnitV)');
    LEGEND_STATES(:,7) = strpad('Pi in component TempRLC (UnitP)');
    LEGEND_STATES(:,8) = strpad('Qo in component TempRLC (UnitQ)');
    LEGEND_CONSTANTS(:,59) = strpad('Rsas in component ParaSys (UnitR)');
    LEGEND_CONSTANTS(:,60) = strpad('Csas in component ParaSys (UnitC)');
    LEGEND_CONSTANTS(:,61) = strpad('Lsas in component ParaSys (UnitL)');
    LEGEND_CONSTANTS(:,62) = strpad('P0sas in component ParaSys (UnitP)');
    LEGEND_CONSTANTS(:,63) = strpad('Q0sas in component ParaSys (UnitQ)');
    LEGEND_ALGEBRAIC(:,33) = strpad('Pi in component TempR (UnitP)');
    LEGEND_STATES(:,9) = strpad('Qo in component TempRLC (UnitQ)');
    LEGEND_CONSTANTS(:,64) = strpad('Rsat in component ParaSys (UnitR)');
    LEGEND_CONSTANTS(:,65) = strpad('Csat in component ParaSys (UnitC)');
    LEGEND_CONSTANTS(:,66) = strpad('Lsat in component ParaSys (UnitL)');
    LEGEND_CONSTANTS(:,67) = strpad('P0sat in component ParaSys (UnitP)');
    LEGEND_CONSTANTS(:,68) = strpad('Q0sat in component ParaSys (UnitQ)');
    LEGEND_ALGEBRAIC(:,29) = strpad('Pi in component TempR (UnitP)');
    LEGEND_ALGEBRAIC(:,26) = strpad('Qo in component TempR (UnitQ)');
    LEGEND_CONSTANTS(:,69) = strpad('Rsar in component ParaSys (UnitR)');
    LEGEND_STATES(:,10) = strpad('Pi in component TempRC (UnitP)');
    LEGEND_ALGEBRAIC(:,28) = strpad('Qo in component TempR (UnitQ)');
    LEGEND_CONSTANTS(:,70) = strpad('Rscp in component ParaSys (UnitR)');
    LEGEND_CONSTANTS(:,71) = strpad('Rsvn in component ParaSys (UnitR)');
    LEGEND_CONSTANTS(:,72) = strpad('Csvn in component ParaSys (UnitC)');
    LEGEND_CONSTANTS(:,73) = strpad('P0svn in component ParaSys (UnitP)');
    LEGEND_STATES(:,11) = strpad('Pi in component TempRLC (UnitP)');
    LEGEND_STATES(:,12) = strpad('Qo in component TempRLC (UnitQ)');
    LEGEND_CONSTANTS(:,74) = strpad('Rpas in component ParaPul (UnitR)');
    LEGEND_CONSTANTS(:,75) = strpad('Cpas in component ParaPul (UnitC)');
    LEGEND_CONSTANTS(:,76) = strpad('Lpas in component ParaPul (UnitL)');
    LEGEND_CONSTANTS(:,77) = strpad('P0pas in component ParaPul (UnitP)');
    LEGEND_CONSTANTS(:,78) = strpad('Q0pas in component ParaPul (UnitQ)');
    LEGEND_ALGEBRAIC(:,34) = strpad('Pi in component TempR (UnitP)');
    LEGEND_STATES(:,13) = strpad('Qo in component TempRLC (UnitQ)');
    LEGEND_CONSTANTS(:,79) = strpad('Rpat in component ParaPul (UnitR)');
    LEGEND_CONSTANTS(:,80) = strpad('Cpat in component ParaPul (UnitC)');
    LEGEND_CONSTANTS(:,81) = strpad('Lpat in component ParaPul (UnitL)');
    LEGEND_CONSTANTS(:,82) = strpad('P0pat in component ParaPul (UnitP)');
    LEGEND_CONSTANTS(:,83) = strpad('Q0pat in component ParaPul (UnitQ)');
    LEGEND_ALGEBRAIC(:,30) = strpad('Pi in component TempR (UnitP)');
    LEGEND_ALGEBRAIC(:,25) = strpad('Qo in component TempR (UnitQ)');
    LEGEND_CONSTANTS(:,84) = strpad('Rpar in component ParaPul (UnitR)');
    LEGEND_STATES(:,14) = strpad('Pi in component TempRC (UnitP)');
    LEGEND_ALGEBRAIC(:,27) = strpad('Qo in component TempR (UnitQ)');
    LEGEND_CONSTANTS(:,85) = strpad('Rpcp in component ParaPul (UnitR)');
    LEGEND_CONSTANTS(:,86) = strpad('Rpvn in component ParaPul (UnitR)');
    LEGEND_CONSTANTS(:,87) = strpad('Cpvn in component ParaPul (UnitC)');
    LEGEND_CONSTANTS(:,88) = strpad('P0pvn in component ParaPul (UnitP)');
    LEGEND_RATES(:,3) = strpad('d/dt V in component TempCDv (UnitV)');
    LEGEND_RATES(:,4) = strpad('d/dt V in component TempCDa (UnitV)');
    LEGEND_RATES(:,5) = strpad('d/dt V in component TempCDv (UnitV)');
    LEGEND_RATES(:,6) = strpad('d/dt V in component TempCDa (UnitV)');
    LEGEND_RATES(:,1) = strpad('d/dt Pi in component TempRLC (UnitP)');
    LEGEND_RATES(:,8) = strpad('d/dt Qo in component TempRLC (UnitQ)');
    LEGEND_RATES(:,7) = strpad('d/dt Pi in component TempRLC (UnitP)');
    LEGEND_RATES(:,9) = strpad('d/dt Qo in component TempRLC (UnitQ)');
    LEGEND_RATES(:,10) = strpad('d/dt Pi in component TempRC (UnitP)');
    LEGEND_RATES(:,2) = strpad('d/dt Pi in component TempRLC (UnitP)');
    LEGEND_RATES(:,12) = strpad('d/dt Qo in component TempRLC (UnitQ)');
    LEGEND_RATES(:,11) = strpad('d/dt Pi in component TempRLC (UnitP)');
    LEGEND_RATES(:,13) = strpad('d/dt Qo in component TempRLC (UnitQ)');
    LEGEND_RATES(:,14) = strpad('d/dt Pi in component TempRC (UnitP)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 350.;
    CONSTANTS(:,2) = 1.0;
    CONSTANTS(:,3) = 5.0;
    CONSTANTS(:,4) = 500;
    CONSTANTS(:,5) = 400.;
    CONSTANTS(:,6) = 1.0;
    CONSTANTS(:,7) = 4.0;
    CONSTANTS(:,8) = 20;
    CONSTANTS(:,9) = 2.5;
    CONSTANTS(:,10) = 0.1;
    CONSTANTS(:,11) = 1.0;
    CONSTANTS(:,12) = 0.3;
    CONSTANTS(:,13) = 0.45;
    CONSTANTS(:,14) = 0.25;
    CONSTANTS(:,15) = 0.15;
    CONSTANTS(:,16) = 0.92;
    CONSTANTS(:,17) = 0.09;
    CONSTANTS(:,18) = 0.25;
    CONSTANTS(:,19) = 0.15;
    CONSTANTS(:,20) = 1.0;
    CONSTANTS(:,21) = 4.0;
    CONSTANTS(:,22) = 1.15;
    CONSTANTS(:,23) = 0.1;
    CONSTANTS(:,24) = 1.0;
    CONSTANTS(:,25) = 10.0;
    CONSTANTS(:,26) = 350.;
    CONSTANTS(:,27) = 400.;
    CONSTANTS(:,28) = 20;
    CONSTANTS(:,29) = 500;
    CONSTANTS(:,30) = 350.;
    CONSTANTS(:,31) = 1.0;
    CONSTANTS(:,32) = 10.0;
    CONSTANTS(:,33) = 500;
    CONSTANTS(:,34) = 400.;
    CONSTANTS(:,35) = 1.0;
    CONSTANTS(:,36) = 4.0;
    CONSTANTS(:,37) = 20;
    CONSTANTS(:,38) = 1.15;
    CONSTANTS(:,39) = 0.1;
    CONSTANTS(:,40) = 1.0;
    CONSTANTS(:,41) = 0.3;
    CONSTANTS(:,42) = 0.45;
    CONSTANTS(:,43) = 0.25;
    CONSTANTS(:,44) = 0.15;
    CONSTANTS(:,45) = 0.92;
    CONSTANTS(:,46) = 0.09;
    CONSTANTS(:,47) = 0.25;
    CONSTANTS(:,48) = 0.15;
    CONSTANTS(:,49) = 1.0;
    CONSTANTS(:,50) = 4.0;
    CONSTANTS(:,51) = 2.5;
    CONSTANTS(:,52) = 0.1;
    CONSTANTS(:,53) = 1.0;
    CONSTANTS(:,54) = 5.0;
    CONSTANTS(:,55) = 350.;
    CONSTANTS(:,56) = 400.;
    CONSTANTS(:,57) = 500;
    CONSTANTS(:,58) = 20;
    CONSTANTS(:,59) = 0.003;
    CONSTANTS(:,60) = 0.08;
    CONSTANTS(:,61) = 0.000062;
    CONSTANTS(:,62) = 100.;
    CONSTANTS(:,63) = 0.;
    CONSTANTS(:,64) = 0.05;
    CONSTANTS(:,65) = 1.6;
    CONSTANTS(:,66) = 0.0017;
    CONSTANTS(:,67) = 100.;
    CONSTANTS(:,68) = 0.;
    CONSTANTS(:,69) = 0.5;
    CONSTANTS(:,70) = 0.52;
    CONSTANTS(:,71) = 0.075;
    CONSTANTS(:,72) = 20.5;
    CONSTANTS(:,73) = 0.;
    CONSTANTS(:,74) = 0.002;
    CONSTANTS(:,75) = 0.18;
    CONSTANTS(:,76) = 0.000052;
    CONSTANTS(:,77) = 30.;
    CONSTANTS(:,78) = 0.;
    CONSTANTS(:,79) = 0.01;
    CONSTANTS(:,80) = 3.8;
    CONSTANTS(:,81) = 0.0017;
    CONSTANTS(:,82) = 30.;
    CONSTANTS(:,83) = 0.;
    CONSTANTS(:,84) = 0.05;
    CONSTANTS(:,85) = 0.25;
    CONSTANTS(:,86) = 0.0006;
    CONSTANTS(:,87) = 20.5;
    CONSTANTS(:,88) = 0.;
    STATES(:,1) = CONSTANTS(:,62);
    STATES(:,2) = CONSTANTS(:,77);
    STATES(:,3) = CONSTANTS(:,4);
    STATES(:,4) = CONSTANTS(:,8);
    STATES(:,5) = CONSTANTS(:,33);
    STATES(:,6) = CONSTANTS(:,37);
    STATES(:,7) = CONSTANTS(:,67);
    STATES(:,8) = CONSTANTS(:,63);
    STATES(:,9) = CONSTANTS(:,68);
    STATES(:,10) = CONSTANTS(:,73);
    STATES(:,11) = CONSTANTS(:,82);
    STATES(:,12) = CONSTANTS(:,78);
    STATES(:,13) = CONSTANTS(:,83);
    STATES(:,14) = CONSTANTS(:,88);
    if (isempty(STATES)), warning('Initial values for states not set');, end
end

function [RATES, ALGEBRAIC] = computeRates(VOI, STATES, CONSTANTS)
    global algebraicVariableCount;
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        ALGEBRAIC = zeros(1, algebraicVariableCount);
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        ALGEBRAIC = zeros(statesRowCount, algebraicVariableCount);
        RATES = zeros(statesRowCount, statesColumnCount);
        utilOnes = ones(statesRowCount, 1);
    end
    RATES(:,8) = ((STATES(:,1) - STATES(:,7)) -  CONSTANTS(:,59).*STATES(:,8))./CONSTANTS(:,61);
    RATES(:,7) = (STATES(:,8) - STATES(:,9))./CONSTANTS(:,65);
    RATES(:,12) = ((STATES(:,2) - STATES(:,11)) -  CONSTANTS(:,74).*STATES(:,12))./CONSTANTS(:,76);
    RATES(:,11) = (STATES(:,12) - STATES(:,13))./CONSTANTS(:,80);
    ALGEBRAIC(:,1) = VOI -  CONSTANTS(:,11).*floor(VOI./CONSTANTS(:,11));
    ALGEBRAIC(:,3) = piecewise({ALGEBRAIC(:,1)>=0.00000&ALGEBRAIC(:,1)<= CONSTANTS(:,12).*CONSTANTS(:,11), 1.00000 - cos(( 3.14159.*ALGEBRAIC(:,1))./( CONSTANTS(:,12).*CONSTANTS(:,11))) , ALGEBRAIC(:,1)> CONSTANTS(:,12).*CONSTANTS(:,11)&ALGEBRAIC(:,1)<= CONSTANTS(:,13).*CONSTANTS(:,11), 1.00000+cos(( 3.14159.*(ALGEBRAIC(:,1) -  CONSTANTS(:,12).*CONSTANTS(:,11)))./( (CONSTANTS(:,13) - CONSTANTS(:,12)).*CONSTANTS(:,11))) , ALGEBRAIC(:,1)> CONSTANTS(:,13).*CONSTANTS(:,11)&ALGEBRAIC(:,1)<CONSTANTS(:,11), 0.00000 }, NaN);
    ALGEBRAIC(:,5) = CONSTANTS(:,10)+( ALGEBRAIC(:,3).*(CONSTANTS(:,9) - CONSTANTS(:,10)))./2.00000;
    ALGEBRAIC(:,7) = CONSTANTS(:,2)+ ALGEBRAIC(:,5).*(STATES(:,3) - CONSTANTS(:,3));
    ALGEBRAIC(:,9) = piecewise({ALGEBRAIC(:,7)>=STATES(:,1), 1.00000 , ALGEBRAIC(:,7)<STATES(:,1), 0.00000 }, NaN);
    ALGEBRAIC(:,11) = piecewise({ALGEBRAIC(:,7)>=STATES(:,1),  CONSTANTS(:,1).*ALGEBRAIC(:,9).*power(abs(ALGEBRAIC(:,7) - STATES(:,1)), 0.500000) , ALGEBRAIC(:,7)<STATES(:,1),  -1.00000.*CONSTANTS(:,1).*ALGEBRAIC(:,9).*power(abs(STATES(:,1) - ALGEBRAIC(:,7)), 0.500000) }, NaN);
    RATES(:,1) = (ALGEBRAIC(:,11) - STATES(:,8))./CONSTANTS(:,60);
    ALGEBRAIC(:,2) = VOI -  CONSTANTS(:,40).*floor(VOI./CONSTANTS(:,40));
    ALGEBRAIC(:,4) = piecewise({ALGEBRAIC(:,2)>=0.00000&ALGEBRAIC(:,2)<= CONSTANTS(:,41).*CONSTANTS(:,40), 1.00000 - cos(( 3.14159.*ALGEBRAIC(:,2))./( CONSTANTS(:,41).*CONSTANTS(:,40))) , ALGEBRAIC(:,2)> CONSTANTS(:,41).*CONSTANTS(:,40)&ALGEBRAIC(:,2)<= CONSTANTS(:,42).*CONSTANTS(:,40), 1.00000+cos(( 3.14159.*(ALGEBRAIC(:,2) -  CONSTANTS(:,41).*CONSTANTS(:,40)))./( (CONSTANTS(:,42) - CONSTANTS(:,41)).*CONSTANTS(:,40))) , ALGEBRAIC(:,2)> CONSTANTS(:,42).*CONSTANTS(:,40)&ALGEBRAIC(:,2)<CONSTANTS(:,40), 0.00000 }, NaN);
    ALGEBRAIC(:,6) = CONSTANTS(:,39)+( ALGEBRAIC(:,4).*(CONSTANTS(:,38) - CONSTANTS(:,39)))./2.00000;
    ALGEBRAIC(:,8) = CONSTANTS(:,31)+ ALGEBRAIC(:,6).*(STATES(:,5) - CONSTANTS(:,32));
    ALGEBRAIC(:,10) = piecewise({ALGEBRAIC(:,8)>=STATES(:,2), 1.00000 , ALGEBRAIC(:,8)<STATES(:,2), 0.00000 }, NaN);
    ALGEBRAIC(:,12) = piecewise({ALGEBRAIC(:,8)>=STATES(:,2),  CONSTANTS(:,30).*ALGEBRAIC(:,10).*power(abs(ALGEBRAIC(:,8) - STATES(:,2)), 0.500000) , ALGEBRAIC(:,8)<STATES(:,2),  -1.00000.*CONSTANTS(:,30).*ALGEBRAIC(:,10).*power(abs(STATES(:,2) - ALGEBRAIC(:,8)), 0.500000) }, NaN);
    RATES(:,2) = (ALGEBRAIC(:,12) - STATES(:,12))./CONSTANTS(:,75);
    ALGEBRAIC(:,13) = VOI -  CONSTANTS(:,11).*floor(VOI./CONSTANTS(:,11));
    ALGEBRAIC(:,15) = piecewise({ALGEBRAIC(:,13)>=0.00000&ALGEBRAIC(:,13)<= (CONSTANTS(:,16)+CONSTANTS(:,17)).*CONSTANTS(:,11) - CONSTANTS(:,11), 1.00000 - cos(( 2.00000.*3.14159.*((ALGEBRAIC(:,13) -  CONSTANTS(:,16).*CONSTANTS(:,11))+CONSTANTS(:,11)))./( CONSTANTS(:,17).*CONSTANTS(:,11))) , ALGEBRAIC(:,13)> (CONSTANTS(:,16)+CONSTANTS(:,17)).*CONSTANTS(:,11) - CONSTANTS(:,11)&ALGEBRAIC(:,13)<= CONSTANTS(:,16).*CONSTANTS(:,11), 0.00000 , ALGEBRAIC(:,13)> CONSTANTS(:,16).*CONSTANTS(:,11)&ALGEBRAIC(:,13)<=CONSTANTS(:,11), 1.00000 - cos(( 2.00000.*3.14159.*(ALGEBRAIC(:,13) -  CONSTANTS(:,16).*CONSTANTS(:,11)))./( CONSTANTS(:,17).*CONSTANTS(:,11))) }, NaN);
    ALGEBRAIC(:,17) = CONSTANTS(:,15)+( ALGEBRAIC(:,15).*(CONSTANTS(:,14) - CONSTANTS(:,15)))./2.00000;
    ALGEBRAIC(:,19) = CONSTANTS(:,6)+ ALGEBRAIC(:,17).*(STATES(:,4) - CONSTANTS(:,7));
    ALGEBRAIC(:,21) = piecewise({ALGEBRAIC(:,19)>=ALGEBRAIC(:,7), 1.00000 , ALGEBRAIC(:,19)<ALGEBRAIC(:,7), 0.00000 }, NaN);
    ALGEBRAIC(:,23) = piecewise({ALGEBRAIC(:,19)>=ALGEBRAIC(:,7),  CONSTANTS(:,5).*ALGEBRAIC(:,21).*power(abs(ALGEBRAIC(:,19) - ALGEBRAIC(:,7)), 0.500000) , ALGEBRAIC(:,19)<ALGEBRAIC(:,7),  -1.00000.*CONSTANTS(:,5).*ALGEBRAIC(:,21).*power(abs(ALGEBRAIC(:,7) - ALGEBRAIC(:,19)), 0.500000) }, NaN);
    RATES(:,3) = ALGEBRAIC(:,23) - ALGEBRAIC(:,11);
    ALGEBRAIC(:,14) = VOI -  CONSTANTS(:,40).*floor(VOI./CONSTANTS(:,40));
    ALGEBRAIC(:,16) = piecewise({ALGEBRAIC(:,14)>=0.00000&ALGEBRAIC(:,14)<= (CONSTANTS(:,45)+CONSTANTS(:,46)).*CONSTANTS(:,40) - CONSTANTS(:,40), 1.00000 - cos(( 2.00000.*3.14159.*((ALGEBRAIC(:,14) -  CONSTANTS(:,45).*CONSTANTS(:,40))+CONSTANTS(:,40)))./( CONSTANTS(:,46).*CONSTANTS(:,40))) , ALGEBRAIC(:,14)> (CONSTANTS(:,45)+CONSTANTS(:,46)).*CONSTANTS(:,40) - CONSTANTS(:,40)&ALGEBRAIC(:,14)<= CONSTANTS(:,45).*CONSTANTS(:,40), 0.00000 , ALGEBRAIC(:,14)> CONSTANTS(:,45).*CONSTANTS(:,40)&ALGEBRAIC(:,14)<=CONSTANTS(:,40), 1.00000 - cos(( 2.00000.*3.14159.*(ALGEBRAIC(:,14) -  CONSTANTS(:,45).*CONSTANTS(:,40)))./( CONSTANTS(:,46).*CONSTANTS(:,40))) }, NaN);
    ALGEBRAIC(:,18) = CONSTANTS(:,44)+( ALGEBRAIC(:,16).*(CONSTANTS(:,43) - CONSTANTS(:,44)))./2.00000;
    ALGEBRAIC(:,20) = CONSTANTS(:,35)+ ALGEBRAIC(:,18).*(STATES(:,6) - CONSTANTS(:,36));
    ALGEBRAIC(:,22) = piecewise({ALGEBRAIC(:,20)>=ALGEBRAIC(:,8), 1.00000 , ALGEBRAIC(:,20)<ALGEBRAIC(:,8), 0.00000 }, NaN);
    ALGEBRAIC(:,24) = piecewise({ALGEBRAIC(:,20)>=ALGEBRAIC(:,8),  CONSTANTS(:,34).*ALGEBRAIC(:,22).*power(abs(ALGEBRAIC(:,20) - ALGEBRAIC(:,8)), 0.500000) , ALGEBRAIC(:,20)<ALGEBRAIC(:,8),  -1.00000.*CONSTANTS(:,34).*ALGEBRAIC(:,22).*power(abs(ALGEBRAIC(:,8) - ALGEBRAIC(:,20)), 0.500000) }, NaN);
    RATES(:,5) = ALGEBRAIC(:,24) - ALGEBRAIC(:,12);
    ALGEBRAIC(:,31) = (STATES(:,14) - ALGEBRAIC(:,19))./CONSTANTS(:,86);
    RATES(:,4) = ALGEBRAIC(:,31) - ALGEBRAIC(:,23);
    ALGEBRAIC(:,32) = (STATES(:,10) - ALGEBRAIC(:,20))./CONSTANTS(:,71);
    RATES(:,6) = ALGEBRAIC(:,32) - ALGEBRAIC(:,24);
    ALGEBRAIC(:,26) = STATES(:,9);
    ALGEBRAIC(:,29) = STATES(:,10)+ CONSTANTS(:,70).*ALGEBRAIC(:,26);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,29)+ CONSTANTS(:,69).*STATES(:,9);
    RATES(:,9) = ((STATES(:,7) - ALGEBRAIC(:,33)) -  CONSTANTS(:,64).*STATES(:,9))./CONSTANTS(:,66);
    ALGEBRAIC(:,28) = ALGEBRAIC(:,26);
    RATES(:,10) = (ALGEBRAIC(:,28) - ALGEBRAIC(:,32))./CONSTANTS(:,72);
    ALGEBRAIC(:,25) = STATES(:,13);
    ALGEBRAIC(:,30) = STATES(:,14)+ CONSTANTS(:,85).*ALGEBRAIC(:,25);
    ALGEBRAIC(:,34) = ALGEBRAIC(:,30)+ CONSTANTS(:,84).*STATES(:,13);
    RATES(:,13) = ((STATES(:,11) - ALGEBRAIC(:,34)) -  CONSTANTS(:,79).*STATES(:,13))./CONSTANTS(:,81);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,25);
    RATES(:,14) = (ALGEBRAIC(:,27) - ALGEBRAIC(:,31))./CONSTANTS(:,87);
   RATES = RATES';
end

% Calculate algebraic variables
function ALGEBRAIC = computeAlgebraic(ALGEBRAIC, CONSTANTS, STATES, VOI)
    statesSize = size(STATES);
    statesColumnCount = statesSize(2);
    if ( statesColumnCount == 1)
        STATES = STATES';
        utilOnes = 1;
    else
        statesRowCount = statesSize(1);
        utilOnes = ones(statesRowCount, 1);
    end
    ALGEBRAIC(:,1) = VOI -  CONSTANTS(:,11).*floor(VOI./CONSTANTS(:,11));
    ALGEBRAIC(:,3) = piecewise({ALGEBRAIC(:,1)>=0.00000&ALGEBRAIC(:,1)<= CONSTANTS(:,12).*CONSTANTS(:,11), 1.00000 - cos(( 3.14159.*ALGEBRAIC(:,1))./( CONSTANTS(:,12).*CONSTANTS(:,11))) , ALGEBRAIC(:,1)> CONSTANTS(:,12).*CONSTANTS(:,11)&ALGEBRAIC(:,1)<= CONSTANTS(:,13).*CONSTANTS(:,11), 1.00000+cos(( 3.14159.*(ALGEBRAIC(:,1) -  CONSTANTS(:,12).*CONSTANTS(:,11)))./( (CONSTANTS(:,13) - CONSTANTS(:,12)).*CONSTANTS(:,11))) , ALGEBRAIC(:,1)> CONSTANTS(:,13).*CONSTANTS(:,11)&ALGEBRAIC(:,1)<CONSTANTS(:,11), 0.00000 }, NaN);
    ALGEBRAIC(:,5) = CONSTANTS(:,10)+( ALGEBRAIC(:,3).*(CONSTANTS(:,9) - CONSTANTS(:,10)))./2.00000;
    ALGEBRAIC(:,7) = CONSTANTS(:,2)+ ALGEBRAIC(:,5).*(STATES(:,3) - CONSTANTS(:,3));
    ALGEBRAIC(:,9) = piecewise({ALGEBRAIC(:,7)>=STATES(:,1), 1.00000 , ALGEBRAIC(:,7)<STATES(:,1), 0.00000 }, NaN);
    ALGEBRAIC(:,11) = piecewise({ALGEBRAIC(:,7)>=STATES(:,1),  CONSTANTS(:,1).*ALGEBRAIC(:,9).*power(abs(ALGEBRAIC(:,7) - STATES(:,1)), 0.500000) , ALGEBRAIC(:,7)<STATES(:,1),  -1.00000.*CONSTANTS(:,1).*ALGEBRAIC(:,9).*power(abs(STATES(:,1) - ALGEBRAIC(:,7)), 0.500000) }, NaN);
    ALGEBRAIC(:,2) = VOI -  CONSTANTS(:,40).*floor(VOI./CONSTANTS(:,40));
    ALGEBRAIC(:,4) = piecewise({ALGEBRAIC(:,2)>=0.00000&ALGEBRAIC(:,2)<= CONSTANTS(:,41).*CONSTANTS(:,40), 1.00000 - cos(( 3.14159.*ALGEBRAIC(:,2))./( CONSTANTS(:,41).*CONSTANTS(:,40))) , ALGEBRAIC(:,2)> CONSTANTS(:,41).*CONSTANTS(:,40)&ALGEBRAIC(:,2)<= CONSTANTS(:,42).*CONSTANTS(:,40), 1.00000+cos(( 3.14159.*(ALGEBRAIC(:,2) -  CONSTANTS(:,41).*CONSTANTS(:,40)))./( (CONSTANTS(:,42) - CONSTANTS(:,41)).*CONSTANTS(:,40))) , ALGEBRAIC(:,2)> CONSTANTS(:,42).*CONSTANTS(:,40)&ALGEBRAIC(:,2)<CONSTANTS(:,40), 0.00000 }, NaN);
    ALGEBRAIC(:,6) = CONSTANTS(:,39)+( ALGEBRAIC(:,4).*(CONSTANTS(:,38) - CONSTANTS(:,39)))./2.00000;
    ALGEBRAIC(:,8) = CONSTANTS(:,31)+ ALGEBRAIC(:,6).*(STATES(:,5) - CONSTANTS(:,32));
    ALGEBRAIC(:,10) = piecewise({ALGEBRAIC(:,8)>=STATES(:,2), 1.00000 , ALGEBRAIC(:,8)<STATES(:,2), 0.00000 }, NaN);
    ALGEBRAIC(:,12) = piecewise({ALGEBRAIC(:,8)>=STATES(:,2),  CONSTANTS(:,30).*ALGEBRAIC(:,10).*power(abs(ALGEBRAIC(:,8) - STATES(:,2)), 0.500000) , ALGEBRAIC(:,8)<STATES(:,2),  -1.00000.*CONSTANTS(:,30).*ALGEBRAIC(:,10).*power(abs(STATES(:,2) - ALGEBRAIC(:,8)), 0.500000) }, NaN);
    ALGEBRAIC(:,13) = VOI -  CONSTANTS(:,11).*floor(VOI./CONSTANTS(:,11));
    ALGEBRAIC(:,15) = piecewise({ALGEBRAIC(:,13)>=0.00000&ALGEBRAIC(:,13)<= (CONSTANTS(:,16)+CONSTANTS(:,17)).*CONSTANTS(:,11) - CONSTANTS(:,11), 1.00000 - cos(( 2.00000.*3.14159.*((ALGEBRAIC(:,13) -  CONSTANTS(:,16).*CONSTANTS(:,11))+CONSTANTS(:,11)))./( CONSTANTS(:,17).*CONSTANTS(:,11))) , ALGEBRAIC(:,13)> (CONSTANTS(:,16)+CONSTANTS(:,17)).*CONSTANTS(:,11) - CONSTANTS(:,11)&ALGEBRAIC(:,13)<= CONSTANTS(:,16).*CONSTANTS(:,11), 0.00000 , ALGEBRAIC(:,13)> CONSTANTS(:,16).*CONSTANTS(:,11)&ALGEBRAIC(:,13)<=CONSTANTS(:,11), 1.00000 - cos(( 2.00000.*3.14159.*(ALGEBRAIC(:,13) -  CONSTANTS(:,16).*CONSTANTS(:,11)))./( CONSTANTS(:,17).*CONSTANTS(:,11))) }, NaN);
    ALGEBRAIC(:,17) = CONSTANTS(:,15)+( ALGEBRAIC(:,15).*(CONSTANTS(:,14) - CONSTANTS(:,15)))./2.00000;
    ALGEBRAIC(:,19) = CONSTANTS(:,6)+ ALGEBRAIC(:,17).*(STATES(:,4) - CONSTANTS(:,7));
    ALGEBRAIC(:,21) = piecewise({ALGEBRAIC(:,19)>=ALGEBRAIC(:,7), 1.00000 , ALGEBRAIC(:,19)<ALGEBRAIC(:,7), 0.00000 }, NaN);
    ALGEBRAIC(:,23) = piecewise({ALGEBRAIC(:,19)>=ALGEBRAIC(:,7),  CONSTANTS(:,5).*ALGEBRAIC(:,21).*power(abs(ALGEBRAIC(:,19) - ALGEBRAIC(:,7)), 0.500000) , ALGEBRAIC(:,19)<ALGEBRAIC(:,7),  -1.00000.*CONSTANTS(:,5).*ALGEBRAIC(:,21).*power(abs(ALGEBRAIC(:,7) - ALGEBRAIC(:,19)), 0.500000) }, NaN);
    ALGEBRAIC(:,14) = VOI -  CONSTANTS(:,40).*floor(VOI./CONSTANTS(:,40));
    ALGEBRAIC(:,16) = piecewise({ALGEBRAIC(:,14)>=0.00000&ALGEBRAIC(:,14)<= (CONSTANTS(:,45)+CONSTANTS(:,46)).*CONSTANTS(:,40) - CONSTANTS(:,40), 1.00000 - cos(( 2.00000.*3.14159.*((ALGEBRAIC(:,14) -  CONSTANTS(:,45).*CONSTANTS(:,40))+CONSTANTS(:,40)))./( CONSTANTS(:,46).*CONSTANTS(:,40))) , ALGEBRAIC(:,14)> (CONSTANTS(:,45)+CONSTANTS(:,46)).*CONSTANTS(:,40) - CONSTANTS(:,40)&ALGEBRAIC(:,14)<= CONSTANTS(:,45).*CONSTANTS(:,40), 0.00000 , ALGEBRAIC(:,14)> CONSTANTS(:,45).*CONSTANTS(:,40)&ALGEBRAIC(:,14)<=CONSTANTS(:,40), 1.00000 - cos(( 2.00000.*3.14159.*(ALGEBRAIC(:,14) -  CONSTANTS(:,45).*CONSTANTS(:,40)))./( CONSTANTS(:,46).*CONSTANTS(:,40))) }, NaN);
    ALGEBRAIC(:,18) = CONSTANTS(:,44)+( ALGEBRAIC(:,16).*(CONSTANTS(:,43) - CONSTANTS(:,44)))./2.00000;
    ALGEBRAIC(:,20) = CONSTANTS(:,35)+ ALGEBRAIC(:,18).*(STATES(:,6) - CONSTANTS(:,36));
    ALGEBRAIC(:,22) = piecewise({ALGEBRAIC(:,20)>=ALGEBRAIC(:,8), 1.00000 , ALGEBRAIC(:,20)<ALGEBRAIC(:,8), 0.00000 }, NaN);
    ALGEBRAIC(:,24) = piecewise({ALGEBRAIC(:,20)>=ALGEBRAIC(:,8),  CONSTANTS(:,34).*ALGEBRAIC(:,22).*power(abs(ALGEBRAIC(:,20) - ALGEBRAIC(:,8)), 0.500000) , ALGEBRAIC(:,20)<ALGEBRAIC(:,8),  -1.00000.*CONSTANTS(:,34).*ALGEBRAIC(:,22).*power(abs(ALGEBRAIC(:,8) - ALGEBRAIC(:,20)), 0.500000) }, NaN);
    ALGEBRAIC(:,31) = (STATES(:,14) - ALGEBRAIC(:,19))./CONSTANTS(:,86);
    ALGEBRAIC(:,32) = (STATES(:,10) - ALGEBRAIC(:,20))./CONSTANTS(:,71);
    ALGEBRAIC(:,26) = STATES(:,9);
    ALGEBRAIC(:,29) = STATES(:,10)+ CONSTANTS(:,70).*ALGEBRAIC(:,26);
    ALGEBRAIC(:,33) = ALGEBRAIC(:,29)+ CONSTANTS(:,69).*STATES(:,9);
    ALGEBRAIC(:,28) = ALGEBRAIC(:,26);
    ALGEBRAIC(:,25) = STATES(:,13);
    ALGEBRAIC(:,30) = STATES(:,14)+ CONSTANTS(:,85).*ALGEBRAIC(:,25);
    ALGEBRAIC(:,34) = ALGEBRAIC(:,30)+ CONSTANTS(:,84).*STATES(:,13);
    ALGEBRAIC(:,27) = ALGEBRAIC(:,25);
end

% Compute result of a piecewise function
function x = piecewise(cases, default)
    set = [0];
    for i = 1:2:length(cases)
        if (length(cases{i+1}) == 1)
            x(cases{i} & ~set,:) = cases{i+1};
        else
            x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
        end
        set = set | cases{i};
        if(set), break, end
    end
    if (length(default) == 1)
        x(~set,:) = default;
    else
        x(~set,:) = default(~set);
    end
end

% Pad out or shorten strings to a set length
function strout = strpad(strin)
    req_length = 160;
    insize = size(strin,2);
    if insize > req_length
        strout = strin(1:req_length);
    else
        strout = [strin, blanks(req_length - insize)];
    end
end
