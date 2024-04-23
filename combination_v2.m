clearvars;
close all;

t0 = 0;
t1 = 48;%hours

% sensitive bacteria, infected bacteria, phage concentration, F = 1 (don't change this)
%initials = [10^6.5, 0, 0, 1];
initials = [10^6.5, 0, 10^5.5, 1];

[t, y] = ode15s(@system, [t0, t1], initials);

total_population = y(:, 1) + y(:, 2);

figure(1);
semilogy(t, y, t, total_population);
ylim([1, 10^7]);
legend(["E: sensitive bacteria", "I: infected bacteria", "C: phage concentration", "F", "E + I: Total bacteria"]);

function dydt = system(t_system, y)
    alpha_gomez = 1.63; %viral particle release rate
    a_gomez = 0.4428; %half-saturation phage density
    beta_gomez = 3.81; %infection rate
    v_gomez = 1.002; %mortality rate
    lambda_gomez = 1.032 * 10^(-2); %viral particle decay rate

    %% copyed from secret_email_correspondence.m

    lambda = 0.11; %bacterial inactivation rate
    gamma = 0.001; %bacterial activation rate
    
    %mu parameters
    alpha = 1.265 * 10^(-4); %growth rate temperature dependence
    T_0 = 2.85; %minimum temperature for growth
    
    %k parameters
    a = 0.11; %decay rate dependence on dew point depression
    L = @(t) 400; %light intensity
    b = 100; %light intensity where the effect on decay is at 50%
    
    %dew point depression threshold
    T_DPD = 2;

    %% the actual model

    mu = @(t) alpha * (T(t) - T_0)^2 * y(4) * (1 - chi(d(t), T_DPD));
    k = @(t) a * (1 + L(t)/(b + L(t))) * d(t) * y(4) * chi(d(t), T_DPD);

    dydt = zeros(4, 1);
    
    %dy/dt
    %If infected bacteria can't reproduce:
    dydt(1) = (mu(t_system) - k(t_system)) * y(1) - beta_gomez * y(1) * y(3) / (a_gomez + y(3));
    %If infected bacteria can reproduce:
    %dydt(1) = (mu(t_system) - k(t_system)) * (y(1) + y(2)) - beta_gomez * y(1) * y(3) / (a_gomez + y(3));
    
    %dI/dt
    dydt(2) = beta_gomez * y(1) * y(3) / (a_gomez + y(3)) - v_gomez * y(2);

    %dC/dt
    dydt(3) = alpha_gomez * y(2) - lambda_gomez * y(3);
    
    %dF/dt
    dydt(4) = -lambda * y(4) * chi(d(t_system), T_DPD) + gamma * (1 - y(4)) * (1 - chi(d(t_system), T_DPD));
end


function temp = T(t)
    %kinda matches the data from fig. 5 and 6

    average_temp = 15;
    temp_variation = 4;
    inoculation_time_offset = 5;

    temp = average_temp + temp_variation .* sin((t - inoculation_time_offset) .* pi/12);
end

function dpd = d(time)
    %{
    relative_humidity = @(t) 50; %going by the dew point calculation, this isn't supposed to be a percentage
    %dew point calculation from supplemental material
    A_1 = 17.625;
    B_1 = 243.04;
    dew_point = @(T, rh) B_1 .* (log(rh ./ 100) + A_1 .* T ./ (B_1 + T)) ./ (A_1 - log(rh ./ 100) - A_1 .* T ./ (B_1 + T)); 
    %dew point depression
    dpd = T(time) - dew_point(T(time), relative_humidity(time));
    %}

    %kinda matches the data from fig. 5 and 6
    dpd = T(time) - 13;
end

function x = chi(dew_point, dew_point_depression_threshold)
    x = zeros(1, length(dew_point));

    for i = 1:length(dew_point)
        if dew_point(i) > dew_point_depression_threshold
            x(i) = 1;
        else
            x(i) = 0;
        end
    end
end