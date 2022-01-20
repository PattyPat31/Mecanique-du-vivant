function [VOI, STATES, ALGEBRAIC, CONSTANTS] = Heartrate()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =0;
end
% There are a total of 0 entries in each of the rate and state variable arrays.
% There are a total of 10 entries in the constant variable array.
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
    LEGEND_CONSTANTS(:,1) = strpad('QLO in component heart_rate_and_stroke_volume (L_per_minute)');
    LEGEND_CONSTANTS(:,2) = strpad('AUR in component heart_rate_and_stroke_volume (dimensionless)');
    LEGEND_CONSTANTS(:,3) = strpad('PRA in component heart_rate_and_stroke_volume (mmHg)');
    LEGEND_CONSTANTS(:,4) = strpad('HMD in component heart_rate_and_stroke_volume (dimensionless)');
    LEGEND_CONSTANTS(:,6) = strpad('AUHR in component effect_of_autonomic_stimulation_on_HR (beats_per_minute)');
    LEGEND_CONSTANTS(:,7) = strpad('PRHR in component effect_of_PRA_on_HR (beats_per_minute)');
    LEGEND_CONSTANTS(:,5) = strpad('PR1LL in component parameter_values (mmHg)');
    LEGEND_CONSTANTS(:,8) = strpad('HDHR in component effect_of_heart_deterioration_on_HR (dimensionless)');
    LEGEND_CONSTANTS(:,9) = strpad('HR in component heart_rate (beats_per_minute)');
    LEGEND_CONSTANTS(:,10) = strpad('SVO in component stroke_volume_output (litre)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 4.9943;
    CONSTANTS(:,2) = 1.30;
    CONSTANTS(:,3) = 0.00852183;
    CONSTANTS(:,4) = 1.0;
    CONSTANTS(:,5) = 0;
    CONSTANTS(:,6) =  72.0000.*CONSTANTS(:,2);
    CONSTANTS(:,7) =  power(CONSTANTS(:,5), 0.500000).*5.00000;
    CONSTANTS(:,8) =  (CONSTANTS(:,4) - 1.00000).*0.500000+1.00000;
    CONSTANTS(:,9) =  (CONSTANTS(:,6)+CONSTANTS(:,7)).*CONSTANTS(:,8);
    CONSTANTS(:,10) = CONSTANTS(:,1)./CONSTANTS(:,9);
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
