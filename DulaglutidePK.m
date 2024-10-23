clear
clc
%% Modeling Dulaglutide PK with a two-compartment model
% V1, C1 correspond to peripheral, subcutaneous (SC) compartment
% V2, C2 correspond to central, plasma compartment

%% Initialize parameters, from Geiser et al. (2016)
D = 1.5e6; % 1.5mg dose, converted to ng to match graph from paper
Ka = 0.00769; % absorption rate constant, /h
Q = 20.1; % inter-compartmental clearance rate, mL/h
CL = 59.3; % clearance rate from central compartment, mL/h
    % Test
    % Q = 20.1e2;
    % Ka = 0.00769e2;
    % CL = 59.3e2; 
F1 = 0.47; % fraction absorbed (bioavailability), unitless
V1 = 3750; % peripheral compartment volume, mL
V2 = 2250; % central compartment volume, mL
p = create_parameter_set([Ka, Q, CL, F1, V1, V2]);

%% Simulate 1.5mg/0.5mL dulaglutide, administered once weekly

% Specify dosing regimen (once-weekly)
tspan = [0, 8]; % 8 weeks
interval = 1; % weekly

% Set initial concentrations in compartments
C1 = D/V1; % sc
C2 = 0; % plasma
C0 = [C1, C2];

% Simulate PK curve
[t_pred,C_pred] = dose_solver(@model,tspan,C0,D,interval,p);

    Cp = C_pred(:,1); % Conc. in peripheral compartment
    Cc = C_pred(:,2); % Conc. in central compartment

%% Plot results
figure(1)
clf
hold on

subplot(2,1,1);
plot(t_pred, Cp, 'LineWidth', 1.5) % Peripheral conc.
title('*S.C.* PK for Once-Weekly 1.5mg S.C. Injection')
xlabel('Time (wk)')
ylabel('Dulaglutide Concentration (ng/mL)')

subplot(2,1,2); 
plot(t_pred, Cc, 'LineWidth', 1.5) % Central conc.
title('Central')
xlabel('Time (wk)')
ylabel('Dulaglutide Concentration (ng/mL)')

% if multiple dose amounts
% for i = 1:length(doses)
%     plot(tspan, pred_C{i}(:,1), 'LineWidth', 1.5) % Peripheral concentrations
%     plot(tspan, pred_C{i}(:,2), 'LineWidth', 1.5) % Central concentrations
% end

hold off
% set(gca, 'YScale', 'log')
title('*Plasma* PK for Once-Weekly 1.5mg S.C. Injection')
xlabel('Time (wk)')
ylabel('Dulaglutide Concentration (ng/mL)')
% legend('Peripheral', 'Central', 'Location', 'best')


%% Functions

function p = create_parameter_set(x)
% create_parameter_set: p = [Ka, Q, CL, F1, V1, V2];
% inputs:
    % x = vector of parameter values
% outputs:
    % pset = parameter set (struct)

p = struct();
p.Ka = x(1);
p.Q = x(2);
p.CL = x(3);
p.F1 = x(4);
p.V1 = x(5);
p.V2 = x(6);
end


function dCdt = model(C,p)
% inputs:
%   C = model state vector
%   p = parameter set (struct)
% outputs:
%   dCdt = model derivative vector

% Pre-allocate output vector
dCdt = zeros(size(C));

% Dependent variables
if length(C) >= 2
    C1 = C(1);  % sc
    C2 = C(2);  % plasma
else
    error('C does not have the expected number of elements.');
end

% Unpack parameters
Ka = p.Ka;
Q = p.Q;
CL = p.CL;
F1 = p.F1;
V1 = p.V1;
V2 = p.V2;

% Rate of change equations
dCdt(1) = -Ka*C1 + (Q/V1)*C2 - (Q/V1)*C1; % sc
dCdt(2) = F1*Ka*C1 + (Q/V2)*C1  - (Q/V2)*C2 - (CL/V2)*C2; % plasma

% Ensure dCdt is a column vector (required for ode45)
dCdt = dCdt(:);  % Convert to a column vector
end


function [t_pred,C_pred] = dose_solver(model,tspan,C0,D,interval,p)
% dose_solver: integrate ODE system with scheduled dosing
% inputs:
    % model = function handle
    % tspan = start and end times (vector)
    % C = initial condition vector
    % D = dosage
    % interval = dosing interval
    % p = model parameters (struct)
% outputs:
    % t = model time (vector)
    % C = integrated model output (vector)

% create output vectors (size is unknown *a priori*)
t_out = [];
C_out = [];

% calculate dosing schedule from inputs
t_start = tspan(1);
t_end = tspan(2);
schedule = t_start:interval:t_end;

% integrate system
for ii = 1:length(schedule)-1 % for each week/dosing point
    tspan = linspace(schedule(ii), schedule(ii+1), 168); % hourly time points (168hrs/week)
    % disp('Initial C0:');
    % disp(C0);
    [t_pred, C_pred] = ode45(@(t,C) model(C,p), tspan, C0);

    % concatenate output and strip final time point (to avoid duplication)
    t_out = cat(1, t_out, t_pred((1:end-1),:));
    C_out = cat(1, C_out, C_pred((1:end-1),:));

    % update initial condition
    C0 = C_pred(end,:);

    % re-dose dulaglutide sc injection
    C0(1) = C0(1) + D/p.V1;

    % C_pred
end

% integrate final interval
if schedule(end) < t_end
    disp('Integrating final interval')
    tspan = linspace(schedule(end), t_end); % 100 time points between last dose & end of study
    [t_pred, C_pred] = ode45(@(t,C) model(C,p), tspan, C0);
    % concatenate output
    t_out = cat(1, t_out, t_pred);
    C_out = cat(1, C_out, C_pred);
    % no re-dosage
end

    t_pred = t_out;
    C_pred = C_out;
end