function [VOI, STATES, ALGEBRAIC, CONSTANTS] = mainFunction()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =22;
end
% There are a total of 10 entries in each of the rate and state variable arrays.
% There are a total of 43 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    tspan = [0, 9];

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
    LEGEND_VOI = strpad('time in component environment (second)');
    LEGEND_CONSTANTS(:,1) = strpad('R_mt in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,2) = strpad('R_av in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,3) = strpad('R_tc in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,4) = strpad('R_pv in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,5) = strpad('R_pul in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,6) = strpad('R_sys in component heart_parameters (kPa_second_per_mL)');
    LEGEND_CONSTANTS(:,7) = strpad('L_tc in component heart_parameters (kPa_second2_per_mL)');
    LEGEND_CONSTANTS(:,8) = strpad('L_pv in component heart_parameters (kPa_second2_per_mL)');
    LEGEND_CONSTANTS(:,9) = strpad('L_mt in component heart_parameters (kPa_second2_per_mL)');
    LEGEND_CONSTANTS(:,10) = strpad('L_av in component heart_parameters (kPa_second2_per_mL)');
    LEGEND_CONSTANTS(:,11) = strpad('V_tot in component heart_parameters (mL)');
    LEGEND_CONSTANTS(:,12) = strpad('P_th in component heart_parameters (kPa)');
    LEGEND_ALGEBRAIC(:,2) = strpad('e_t in component driver_function (dimensionless)');
    LEGEND_CONSTANTS(:,13) = strpad('A in component driver_function (dimensionless)');
    LEGEND_CONSTANTS(:,14) = strpad('B in component driver_function (per_second2)');
    LEGEND_CONSTANTS(:,15) = strpad('C in component driver_function (second)');
    LEGEND_ALGEBRAIC(:,1) = strpad('tau in component driver_function (second)');
    LEGEND_CONSTANTS(:,16) = strpad('period in component driver_function (second)');
    LEGEND_ALGEBRAIC(:,3) = strpad('V_pcd in component pericardium (mL)');
    LEGEND_ALGEBRAIC(:,4) = strpad('P_pcd in component pericardium (kPa)');
    LEGEND_ALGEBRAIC(:,5) = strpad('P_peri in component pericardium (kPa)');
    LEGEND_STATES(:,1) = strpad('V_lv in component left_ventricle (mL)');
    LEGEND_STATES(:,2) = strpad('V_rv in component right_ventricle (mL)');
    LEGEND_CONSTANTS(:,17) = strpad('P_0_pcd in component pericardium (kPa)');
    LEGEND_CONSTANTS(:,18) = strpad('V_0_pcd in component pericardium (mL)');
    LEGEND_CONSTANTS(:,19) = strpad('lambda_pcd in component pericardium (per_mL)');
    LEGEND_ALGEBRAIC(:,7) = strpad('V_lvf in component left_ventricle (mL)');
    LEGEND_ALGEBRAIC(:,10) = strpad('P_lvf in component left_ventricle (kPa)');
    LEGEND_ALGEBRAIC(:,11) = strpad('P_lv in component left_ventricle (kPa)');
    LEGEND_ALGEBRAIC(:,6) = strpad('V_spt in component septum (mL)');
    LEGEND_ALGEBRAIC(:,8) = strpad('P_es_lvf in component lvf_calculator (kPa)');
    LEGEND_ALGEBRAIC(:,9) = strpad('P_ed_lvf in component lvf_calculator (kPa)');
    LEGEND_ALGEBRAIC(:,19) = strpad('P_pu in component pulmonary_vein (kPa)');
    LEGEND_ALGEBRAIC(:,18) = strpad('P_ao in component aorta (kPa)');
    LEGEND_CONSTANTS(:,20) = strpad('E_es_lvf in component lvf_calculator (kPa_per_mL)');
    LEGEND_CONSTANTS(:,21) = strpad('lambda_lvf in component lvf_calculator (per_mL)');
    LEGEND_CONSTANTS(:,22) = strpad('P_0_lvf in component lvf_calculator (kPa)');
    LEGEND_STATES(:,3) = strpad('Q_mt in component flow (mL_per_second)');
    LEGEND_STATES(:,4) = strpad('Q_av in component flow (mL_per_second)');
    LEGEND_CONSTANTS(:,23) = strpad('V_d_lvf in component lvf_calculator (mL)');
    LEGEND_CONSTANTS(:,24) = strpad('V_0_lvf in component lvf_calculator (mL)');
    LEGEND_ALGEBRAIC(:,12) = strpad('V_rvf in component right_ventricle (mL)');
    LEGEND_ALGEBRAIC(:,15) = strpad('P_rvf in component right_ventricle (kPa)');
    LEGEND_ALGEBRAIC(:,16) = strpad('P_rv in component right_ventricle (kPa)');
    LEGEND_ALGEBRAIC(:,13) = strpad('P_es_rvf in component rvf_calculator (kPa)');
    LEGEND_ALGEBRAIC(:,14) = strpad('P_ed_rvf in component rvf_calculator (kPa)');
    LEGEND_ALGEBRAIC(:,17) = strpad('P_pa in component pulmonary_artery (kPa)');
    LEGEND_ALGEBRAIC(:,20) = strpad('P_vc in component vena_cava (kPa)');
    LEGEND_CONSTANTS(:,25) = strpad('E_es_rvf in component rvf_calculator (kPa_per_mL)');
    LEGEND_CONSTANTS(:,26) = strpad('lambda_rvf in component rvf_calculator (per_mL)');
    LEGEND_CONSTANTS(:,27) = strpad('P_0_rvf in component rvf_calculator (kPa)');
    LEGEND_STATES(:,5) = strpad('Q_tc in component flow (mL_per_second)');
    LEGEND_STATES(:,6) = strpad('Q_pv in component flow (mL_per_second)');
    LEGEND_CONSTANTS(:,28) = strpad('V_d_rvf in component rvf_calculator (mL)');
    LEGEND_CONSTANTS(:,29) = strpad('V_0_rvf in component rvf_calculator (mL)');
    LEGEND_CONSTANTS(:,30) = strpad('E_es_spt in component septum (kPa_per_mL)');
    LEGEND_CONSTANTS(:,31) = strpad('V_d_spt in component septum (mL)');
    LEGEND_CONSTANTS(:,32) = strpad('P_0_spt in component septum (kPa)');
    LEGEND_CONSTANTS(:,33) = strpad('lambda_spt in component septum (per_mL)');
    LEGEND_CONSTANTS(:,34) = strpad('V_0_spt in component septum (mL)');
    LEGEND_CONSTANTS(:,35) = strpad('one in component septum (dimensionless)');
    LEGEND_CONSTANTS(:,36) = strpad('E_es_pa in component pulmonary_artery (kPa_per_mL)');
    LEGEND_STATES(:,7) = strpad('V_pa in component pulmonary_artery (mL)');
    LEGEND_CONSTANTS(:,37) = strpad('V_d_pa in component pulmonary_artery (mL)');
    LEGEND_ALGEBRAIC(:,21) = strpad('Q_pul in component flow (mL_per_second)');
    LEGEND_CONSTANTS(:,38) = strpad('E_es_pu in component pulmonary_vein (kPa_per_mL)');
    LEGEND_STATES(:,8) = strpad('V_pu in component pulmonary_vein (mL)');
    LEGEND_CONSTANTS(:,39) = strpad('V_d_pu in component pulmonary_vein (mL)');
    LEGEND_CONSTANTS(:,40) = strpad('E_es_ao in component aorta (kPa_per_mL)');
    LEGEND_STATES(:,9) = strpad('V_ao in component aorta (mL)');
    LEGEND_CONSTANTS(:,41) = strpad('V_d_ao in component aorta (mL)');
    LEGEND_ALGEBRAIC(:,22) = strpad('Q_sys in component flow (mL_per_second)');
    LEGEND_CONSTANTS(:,42) = strpad('E_es_vc in component vena_cava (kPa_per_mL)');
    LEGEND_STATES(:,10) = strpad('V_vc in component vena_cava (mL)');
    LEGEND_CONSTANTS(:,43) = strpad('V_d_vc in component vena_cava (mL)');
    LEGEND_RATES(:,1) = strpad('d/dt V_lv in component left_ventricle (mL)');
    LEGEND_RATES(:,2) = strpad('d/dt V_rv in component right_ventricle (mL)');
    LEGEND_RATES(:,7) = strpad('d/dt V_pa in component pulmonary_artery (mL)');
    LEGEND_RATES(:,8) = strpad('d/dt V_pu in component pulmonary_vein (mL)');
    LEGEND_RATES(:,9) = strpad('d/dt V_ao in component aorta (mL)');
    LEGEND_RATES(:,10) = strpad('d/dt V_vc in component vena_cava (mL)');
    LEGEND_RATES(:,3) = strpad('d/dt Q_mt in component flow (mL_per_second)');
    LEGEND_RATES(:,4) = strpad('d/dt Q_av in component flow (mL_per_second)');
    LEGEND_RATES(:,5) = strpad('d/dt Q_tc in component flow (mL_per_second)');
    LEGEND_RATES(:,6) = strpad('d/dt Q_pv in component flow (mL_per_second)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 0.0158;
    CONSTANTS(:,2) = 0.0180;
    CONSTANTS(:,3) = 0.0237;
    CONSTANTS(:,4) = 0.0055;
    CONSTANTS(:,5) = 0.1552;
    CONSTANTS(:,6) = 1.0889;
    CONSTANTS(:,7) = 8.0093e-5;
    CONSTANTS(:,8) = 1.4868e-4;
    CONSTANTS(:,9) = 7.6968e-5;
    CONSTANTS(:,10) = 1.2189e-4;
    CONSTANTS(:,11) = 5.5;
    CONSTANTS(:,12) = -4;
    CONSTANTS(:,13) = 1;
    CONSTANTS(:,14) = 80;
    CONSTANTS(:,15) = 0.375;
    CONSTANTS(:,16) = 0.75;
    STATES(:,1) = 94.6812;
    STATES(:,2) = 90.7302;
    CONSTANTS(:,17) = 0.5003;
    CONSTANTS(:,18) = 200;
    CONSTANTS(:,19) = 0.03;
    CONSTANTS(:,20) = 2.8798;
    CONSTANTS(:,21) = 0.033;
    CONSTANTS(:,22) = 0.1203;
    STATES(:,3) = 245.5813;
    STATES(:,4) = 0;
    CONSTANTS(:,23) = 0;
    CONSTANTS(:,24) = 0;
    CONSTANTS(:,25) = 0.585;
    CONSTANTS(:,26) = 0.023;
    CONSTANTS(:,27) = 0.2157;
    STATES(:,5) = 190.0661;
    STATES(:,6) = 0;
    CONSTANTS(:,28) = 0;
    CONSTANTS(:,29) = 0;
    CONSTANTS(:,30) = 48.754;
    CONSTANTS(:,31) = 2;
    CONSTANTS(:,32) = 1.1101;
    CONSTANTS(:,33) = 0.435;
    CONSTANTS(:,34) = 2;
    CONSTANTS(:,35) = 1;
    CONSTANTS(:,36) = 0.369;
    STATES(:,7) = 43.0123;
    CONSTANTS(:,37) = 0;
    CONSTANTS(:,38) = 0.0073;
    STATES(:,8) = 808.4579;
    CONSTANTS(:,39) = 0;
    CONSTANTS(:,40) = 0.6913;
    STATES(:,9) = 133.3381;
    CONSTANTS(:,41) = 0;
    CONSTANTS(:,42) = 0.0059;
    STATES(:,10) = 329.7803;
    CONSTANTS(:,43) = 0;
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
    RATES(:,1) = piecewise({STATES(:,3)<0.00000&STATES(:,4)<0.00000, 0.00000 , STATES(:,3)<0.00000,  - STATES(:,4) , STATES(:,4)<0.00000, STATES(:,3) }, STATES(:,3) - STATES(:,4));
    RATES(:,2) = piecewise({STATES(:,5)<0.00000&STATES(:,6)<0.00000, 0.00000 , STATES(:,5)<0.00000,  - STATES(:,6) , STATES(:,6)<0.00000, STATES(:,5) }, STATES(:,5) - STATES(:,6));
    ALGEBRAIC(:,3) = STATES(:,1)+STATES(:,2);
    ALGEBRAIC(:,4) =  CONSTANTS(:,17).*(exp( CONSTANTS(:,19).*(ALGEBRAIC(:,3) - CONSTANTS(:,18))) - 1.00000);
    ALGEBRAIC(:,5) = ALGEBRAIC(:,4)+CONSTANTS(:,12);
    ALGEBRAIC(:,1) = piecewise({VOI<=CONSTANTS(:,16), VOI , VOI<= CONSTANTS(:,16).*2.00000, VOI - CONSTANTS(:,16) , VOI<= CONSTANTS(:,16).*3.00000, VOI -  CONSTANTS(:,16).*2.00000 , VOI<= CONSTANTS(:,16).*4.00000, VOI -  CONSTANTS(:,16).*3.00000 , VOI<= CONSTANTS(:,16).*5.00000, VOI -  CONSTANTS(:,16).*4.00000 , VOI<= CONSTANTS(:,16).*6.00000, VOI -  CONSTANTS(:,16).*5.00000 , VOI<= CONSTANTS(:,16).*7.00000, VOI -  CONSTANTS(:,16).*6.00000 , VOI<= CONSTANTS(:,16).*8.00000, VOI -  CONSTANTS(:,16).*7.00000 , VOI<= CONSTANTS(:,16).*9.00000, VOI -  CONSTANTS(:,16).*8.00000 , VOI<= CONSTANTS(:,16).*10.0000, VOI -  CONSTANTS(:,16).*9.00000 , VOI<= CONSTANTS(:,16).*11.0000, VOI -  CONSTANTS(:,16).*10.0000 , VOI<= CONSTANTS(:,16).*12.0000, VOI -  CONSTANTS(:,16).*11.0000 , VOI<= CONSTANTS(:,16).*13.0000, VOI -  CONSTANTS(:,16).*12.0000 }, NaN);
    ALGEBRAIC(:,2) =  CONSTANTS(:,13).*exp(  - CONSTANTS(:,14).*power(ALGEBRAIC(:,1) - CONSTANTS(:,15), 2.00000));
    [CONSTANTS, STATES, ALGEBRAIC] = rootfind_0(VOI, CONSTANTS, STATES, ALGEBRAIC);
    ALGEBRAIC(:,7) = STATES(:,1) - ALGEBRAIC(:,6);
    ALGEBRAIC(:,8) =  CONSTANTS(:,20).*(ALGEBRAIC(:,7) - CONSTANTS(:,23));
    ALGEBRAIC(:,9) =  CONSTANTS(:,22).*(exp( CONSTANTS(:,21).*(ALGEBRAIC(:,7) - CONSTANTS(:,24))) - 1.00000);
    ALGEBRAIC(:,10) =  ALGEBRAIC(:,2).*ALGEBRAIC(:,8)+ (1.00000 - ALGEBRAIC(:,2)).*ALGEBRAIC(:,9);
    ALGEBRAIC(:,11) = ALGEBRAIC(:,10)+ALGEBRAIC(:,5);
    ALGEBRAIC(:,18) =  CONSTANTS(:,40).*(STATES(:,9) - CONSTANTS(:,41));
    RATES(:,4) = piecewise({ALGEBRAIC(:,11) - ALGEBRAIC(:,18)<0.00000&STATES(:,4)<0.00000, 0.00000 }, ((ALGEBRAIC(:,11) - ALGEBRAIC(:,18)) -  STATES(:,4).*CONSTANTS(:,2))./CONSTANTS(:,10));
    ALGEBRAIC(:,12) = STATES(:,2)+ALGEBRAIC(:,6);
    ALGEBRAIC(:,13) =  CONSTANTS(:,25).*(ALGEBRAIC(:,12) - CONSTANTS(:,28));
    ALGEBRAIC(:,14) =  CONSTANTS(:,27).*(exp( CONSTANTS(:,26).*(ALGEBRAIC(:,12) - CONSTANTS(:,29))) - 1.00000);
    ALGEBRAIC(:,15) =  ALGEBRAIC(:,2).*ALGEBRAIC(:,13)+ (1.00000 - ALGEBRAIC(:,2)).*ALGEBRAIC(:,14);
    ALGEBRAIC(:,16) = ALGEBRAIC(:,15)+ALGEBRAIC(:,5);
    ALGEBRAIC(:,17) =  CONSTANTS(:,36).*(STATES(:,7) - CONSTANTS(:,37))+CONSTANTS(:,12);
    RATES(:,6) = piecewise({ALGEBRAIC(:,16) - ALGEBRAIC(:,17)<0.00000&STATES(:,6)<0.00000, 0.00000 }, ((ALGEBRAIC(:,16) - ALGEBRAIC(:,17)) -  STATES(:,6).*CONSTANTS(:,4))./CONSTANTS(:,8));
    ALGEBRAIC(:,19) =  CONSTANTS(:,38).*(STATES(:,8) - CONSTANTS(:,39))+CONSTANTS(:,12);
    RATES(:,3) = piecewise({ALGEBRAIC(:,19) - ALGEBRAIC(:,11)<0.00000&STATES(:,3)<0.00000, 0.00000 }, ((ALGEBRAIC(:,19) - ALGEBRAIC(:,11)) -  STATES(:,3).*CONSTANTS(:,1))./CONSTANTS(:,9));
    ALGEBRAIC(:,20) =  CONSTANTS(:,42).*(STATES(:,10) - CONSTANTS(:,43));
    RATES(:,5) = piecewise({ALGEBRAIC(:,20) - ALGEBRAIC(:,16)<0.00000&STATES(:,5)<0.00000, 0.00000 }, ((ALGEBRAIC(:,20) - ALGEBRAIC(:,16)) -  STATES(:,5).*CONSTANTS(:,3))./CONSTANTS(:,7));
    ALGEBRAIC(:,21) = (ALGEBRAIC(:,17) - ALGEBRAIC(:,19))./CONSTANTS(:,5);
    RATES(:,7) = piecewise({STATES(:,6)<0.00000,  - ALGEBRAIC(:,21) }, STATES(:,6) - ALGEBRAIC(:,21));
    RATES(:,8) = piecewise({STATES(:,3)<0.00000, ALGEBRAIC(:,21) }, ALGEBRAIC(:,21) - STATES(:,3));
    ALGEBRAIC(:,22) = (ALGEBRAIC(:,18) - ALGEBRAIC(:,20))./CONSTANTS(:,6);
    RATES(:,9) = piecewise({STATES(:,4)<0.00000,  - ALGEBRAIC(:,22) }, STATES(:,4) - ALGEBRAIC(:,22));
    RATES(:,10) = piecewise({STATES(:,5)<0.00000, ALGEBRAIC(:,22) }, ALGEBRAIC(:,22) - STATES(:,5));
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
    ALGEBRAIC(:,3) = STATES(:,1)+STATES(:,2);
    ALGEBRAIC(:,4) =  CONSTANTS(:,17).*(exp( CONSTANTS(:,19).*(ALGEBRAIC(:,3) - CONSTANTS(:,18))) - 1.00000);
    ALGEBRAIC(:,5) = ALGEBRAIC(:,4)+CONSTANTS(:,12);
    ALGEBRAIC(:,1) = piecewise({VOI<=CONSTANTS(:,16), VOI , VOI<= CONSTANTS(:,16).*2.00000, VOI - CONSTANTS(:,16) , VOI<= CONSTANTS(:,16).*3.00000, VOI -  CONSTANTS(:,16).*2.00000 , VOI<= CONSTANTS(:,16).*4.00000, VOI -  CONSTANTS(:,16).*3.00000 , VOI<= CONSTANTS(:,16).*5.00000, VOI -  CONSTANTS(:,16).*4.00000 , VOI<= CONSTANTS(:,16).*6.00000, VOI -  CONSTANTS(:,16).*5.00000 , VOI<= CONSTANTS(:,16).*7.00000, VOI -  CONSTANTS(:,16).*6.00000 , VOI<= CONSTANTS(:,16).*8.00000, VOI -  CONSTANTS(:,16).*7.00000 , VOI<= CONSTANTS(:,16).*9.00000, VOI -  CONSTANTS(:,16).*8.00000 , VOI<= CONSTANTS(:,16).*10.0000, VOI -  CONSTANTS(:,16).*9.00000 , VOI<= CONSTANTS(:,16).*11.0000, VOI -  CONSTANTS(:,16).*10.0000 , VOI<= CONSTANTS(:,16).*12.0000, VOI -  CONSTANTS(:,16).*11.0000 , VOI<= CONSTANTS(:,16).*13.0000, VOI -  CONSTANTS(:,16).*12.0000 }, NaN);
    ALGEBRAIC(:,2) =  CONSTANTS(:,13).*exp(  - CONSTANTS(:,14).*power(ALGEBRAIC(:,1) - CONSTANTS(:,15), 2.00000));
    ALGEBRAIC(:,7) = STATES(:,1) - ALGEBRAIC(:,6);
    ALGEBRAIC(:,8) =  CONSTANTS(:,20).*(ALGEBRAIC(:,7) - CONSTANTS(:,23));
    ALGEBRAIC(:,9) =  CONSTANTS(:,22).*(exp( CONSTANTS(:,21).*(ALGEBRAIC(:,7) - CONSTANTS(:,24))) - 1.00000);
    ALGEBRAIC(:,10) =  ALGEBRAIC(:,2).*ALGEBRAIC(:,8)+ (1.00000 - ALGEBRAIC(:,2)).*ALGEBRAIC(:,9);
    ALGEBRAIC(:,11) = ALGEBRAIC(:,10)+ALGEBRAIC(:,5);
    ALGEBRAIC(:,18) =  CONSTANTS(:,40).*(STATES(:,9) - CONSTANTS(:,41));
    ALGEBRAIC(:,12) = STATES(:,2)+ALGEBRAIC(:,6);
    ALGEBRAIC(:,13) =  CONSTANTS(:,25).*(ALGEBRAIC(:,12) - CONSTANTS(:,28));
    ALGEBRAIC(:,14) =  CONSTANTS(:,27).*(exp( CONSTANTS(:,26).*(ALGEBRAIC(:,12) - CONSTANTS(:,29))) - 1.00000);
    ALGEBRAIC(:,15) =  ALGEBRAIC(:,2).*ALGEBRAIC(:,13)+ (1.00000 - ALGEBRAIC(:,2)).*ALGEBRAIC(:,14);
    ALGEBRAIC(:,16) = ALGEBRAIC(:,15)+ALGEBRAIC(:,5);
    ALGEBRAIC(:,17) =  CONSTANTS(:,36).*(STATES(:,7) - CONSTANTS(:,37))+CONSTANTS(:,12);
    ALGEBRAIC(:,19) =  CONSTANTS(:,38).*(STATES(:,8) - CONSTANTS(:,39))+CONSTANTS(:,12);
    ALGEBRAIC(:,20) =  CONSTANTS(:,42).*(STATES(:,10) - CONSTANTS(:,43));
    ALGEBRAIC(:,21) = (ALGEBRAIC(:,17) - ALGEBRAIC(:,19))./CONSTANTS(:,5);
    ALGEBRAIC(:,22) = (ALGEBRAIC(:,18) - ALGEBRAIC(:,20))./CONSTANTS(:,6);
end

% Functions required for solving differential algebraic equation
function [CONSTANTS, STATES, ALGEBRAIC] = rootfind_0(VOI, CONSTANTS_IN, STATES_IN, ALGEBRAIC_IN)
    CONSTANTS = CONSTANTS_IN;
    STATES = STATES_IN;
    ALGEBRAIC = ALGEBRAIC_IN;
    global initialGuess_0;
    if (length(initialGuess_0) ~= 1), initialGuess_0 = 0.1;, end
    options = optimset('Display', 'off', 'TolX', 1E-6);
    if length(VOI) == 1
        residualfn = @(algebraicCandidate)residualSN_0(algebraicCandidate, ALGEBRAIC, VOI, CONSTANTS, STATES);
        ALGEBRAIC(:,6) = fsolve(residualfn, initialGuess_0, options);
        initialGuess_0 = ALGEBRAIC(:,6);
    else
        SET_ALGEBRAIC(:,6) = logical(1);
        for i=1:length(VOI)
            residualfn = @(algebraicCandidate)residualSN_0(algebraicCandidate, ALGEBRAIC(i,:), VOI(i), CONSTANTS, STATES(i,:));
            TEMP_ALGEBRAIC(:,6) = fsolve(residualfn, initialGuess_0, options);
            ALGEBRAIC(i,SET_ALGEBRAIC) = TEMP_ALGEBRAIC(SET_ALGEBRAIC);
            initialGuess_0 = TEMP_ALGEBRAIC(:,6);
        end
    end
end

function resid = residualSN_0(algebraicCandidate, ALGEBRAIC, VOI, CONSTANTS, STATES)
    ALGEBRAIC(:,6) = algebraicCandidate;
    resid = (0.00000) - (((( ALGEBRAIC(:,2).*CONSTANTS(:,30).*(ALGEBRAIC(:,6) - CONSTANTS(:,31))+ (CONSTANTS(:,35) - ALGEBRAIC(:,2)).*CONSTANTS(:,32).*(exp( CONSTANTS(:,33).*(ALGEBRAIC(:,6) - CONSTANTS(:,34))) - CONSTANTS(:,35))) -  ALGEBRAIC(:,2).*CONSTANTS(:,20).*(STATES(:,1) - ALGEBRAIC(:,6))) -  (1.00000 - ALGEBRAIC(:,2)).*CONSTANTS(:,22).*(exp( CONSTANTS(:,21).*(STATES(:,1) - ALGEBRAIC(:,6))) - 1.00000))+ ALGEBRAIC(:,2).*CONSTANTS(:,25).*(STATES(:,2)+ALGEBRAIC(:,6))+ (1.00000 - ALGEBRAIC(:,2)).*CONSTANTS(:,27).*(exp( CONSTANTS(:,26).*(STATES(:,2)+ALGEBRAIC(:,6))) - 1.00000));
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
