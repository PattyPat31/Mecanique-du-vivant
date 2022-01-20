function [VOI, STATES, ALGEBRAIC, CONSTANTS] = coeur()
    % This is the "main function".  In Matlab, things work best if you rename this function to match the filename.
   [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel();
end

function [algebraicVariableCount] = getAlgebraicVariableCount() 
    % Used later when setting a global variable with the number of algebraic variables.
    % Note: This is not the "main method".  
    algebraicVariableCount =2;
end
% There are a total of 1 entries in each of the rate and state variable arrays.
% There are a total of 7 entries in the constant variable array.
%

function [VOI, STATES, ALGEBRAIC, CONSTANTS] = solveModel()
    % Create ALGEBRAIC of correct size
    global algebraicVariableCount;  algebraicVariableCount = getAlgebraicVariableCount();
    % Initialise constants and state variables
    [INIT_STATES, CONSTANTS] = initConsts;

    % Set timespan to solve over 
    %tspan = [0, 10];
    tspan=linspace(0,10,500);

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
    hold on;
    plot(VOI,ALGEBRAIC);
    xlabel(LEGEND_VOI);
    l = legend(LEGEND_STATES,LEGEND_ALGEBRAIC(1,:),LEGEND_ALGEBRAIC(2,:));
    set(l,'Interpreter','none');
    hold off;
    
end

function [LEGEND_STATES, LEGEND_ALGEBRAIC, LEGEND_VOI, LEGEND_CONSTANTS] = createLegends()
    LEGEND_STATES = ''; LEGEND_ALGEBRAIC = ''; LEGEND_VOI = ''; LEGEND_CONSTANTS = '';
    LEGEND_VOI = strpad('time in component Environment (s)');
    LEGEND_CONSTANTS(:,1) = strpad('HR in component Environment (ratepm)');
    LEGEND_CONSTANTS(:,6) = strpad('hrf in component Environment (Hz)');
    LEGEND_CONSTANTS(:,2) = strpad('PRint in component LVTiming (s)');
    LEGEND_ALGEBRAIC(:,1) = strpad('beattime in component LVTiming (s)');
    LEGEND_CONSTANTS(:,3) = strpad('Esys in component LVElastanceFunction (elastance)');
    LEGEND_CONSTANTS(:,4) = strpad('Edia in component LVElastanceFunction (elastance)');
    LEGEND_CONSTANTS(:,5) = strpad('TsK in component LVElastanceFunction (s)');
    LEGEND_CONSTANTS(:,7) = strpad('Ts in component LVElastanceFunction (s)');
    LEGEND_ALGEBRAIC(:,2) = strpad('E_LV in component LVElastanceFunction (elastance)');
    LEGEND_STATES(:,1) = strpad('dummy in component dummy (dimensionless)');
    LEGEND_RATES(:,1) = strpad('d/dt dummy in component dummy (dimensionless)');
    LEGEND_STATES  = LEGEND_STATES';
    LEGEND_ALGEBRAIC = LEGEND_ALGEBRAIC';
    LEGEND_RATES = LEGEND_RATES';
    LEGEND_CONSTANTS = LEGEND_CONSTANTS';
end

function [STATES, CONSTANTS] = initConsts()
    VOI = 0; CONSTANTS = []; STATES = []; ALGEBRAIC = [];
    CONSTANTS(:,1) = 70;
    CONSTANTS(:,2) = 0.00012;
    CONSTANTS(:,3) = 5.6;
    CONSTANTS(:,4) = 0.19;
    CONSTANTS(:,5) = 0.35;
    STATES(:,1) = 10;
    CONSTANTS(:,6) = CONSTANTS(:,1)./60.0000;
    CONSTANTS(:,7) =  CONSTANTS(:,5).*power(( 1.00000.*CONSTANTS(:,6)), 1.0 ./ 2);
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
    RATES(:,1) =  STATES(:,1).* - 3.00000;
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
    ALGEBRAIC=[0.0*VOI, 0.0*VOI];
    ALGEBRAIC(:,1) = (VOI -  floor(VOI./CONSTANTS(:,6)).*CONSTANTS(:,6)) - CONSTANTS(:,2);
    ALGEBRAIC(:,2) = piecewise({ALGEBRAIC(:,1)>=0.00000&ALGEBRAIC(:,1)<=CONSTANTS(:,7), CONSTANTS(:,4)+( (CONSTANTS(:,3) - CONSTANTS(:,4)).*(1.00000 - cos((  pi.*ALGEBRAIC(:,1))./CONSTANTS(:,7))))./2.00000 , ALGEBRAIC(:,1)< 1.50000.*CONSTANTS(:,7)&ALGEBRAIC(:,1)>=CONSTANTS(:,7), CONSTANTS(:,4)+( (CONSTANTS(:,3) - CONSTANTS(:,4)).*(1.00000+cos(( 2.00000.* pi.*(ALGEBRAIC(:,1) - CONSTANTS(:,7)))./CONSTANTS(:,7))))./2.00000 }, CONSTANTS(:,4));
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