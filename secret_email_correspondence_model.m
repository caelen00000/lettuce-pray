clearvars;
close all;

initial_population = 10^6.5;

initials = [initial_population, 1];

t0 = 0;
t1 = 175;

[t, y] = ode45(@system, [t0, t1], initials);

figure(1);
semilogy(t, y);
legend(["pop", "F"]);

times = t0:0.01:t1;

lambda = 0.11; %bacterial inactivation rate
gamma = 0.001; %bacterial activation rate

%mu parameters
alpha = 1.265 * 10^(-4); %scale factor
T_0 = 2.85; %minimum temperature for growth

%k parameters
a = 0.11; %scale factor
L = @(t) 400; %light intensity
b = 100; %light intensity where the effect on decay is at 50%

%dew point depression threshold
T_DPD = 2;

mu = @(t) (alpha .* (T(t) - T_0).^2) .* (1 - chi(d(t), T_DPD));
k = @(t) a .* (1 + L(t)/(b + L(t))) .* d(t) .* chi(d(t), T_DPD);

figure(2);
plot(times, T(times), times, d(times));

figure(3);
plot(times, mu(times), times, k(times));

function dydt = system(t_system, y)
    lambda = 0.11; %bacterial inactivation rate
    gamma = 0.001; %bacterial activation rate
    
    %mu parameters
    alpha = 1.265 * 10^(-4); %scale factor
    T_0 = 2.85; %minimum temperature for growth
    
    %k parameters
    a = 0.11; %scale factor
    L = @(t) 400; %light intensity
    b = 100; %light intensity where the effect on decay is at 50%
    
    %dew point depression threshold
    T_DPD = 2;

    mu = @(t) alpha * (T(t) - T_0)^2 * y(2) * (1 - chi(d(t), T_DPD));
    k = @(t) a * (1 + L(t)/(b + L(t))) * d(t) * y(2) * chi(d(t), T_DPD);

    dydt = zeros(2, 1);
    dydt(1) = (mu(t_system) - k(t_system)) * y(1); %dy/dt
    dydt(2) = -lambda * y(2) * chi(d(t_system), T_DPD) + gamma * (1 - y(2)) * (1 - chi(d(t_system), T_DPD));%dF/dt
end

function temp = T(t)
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