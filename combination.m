clearvars;
close all;

t0 = 0;
t1 = 4;

initials = [0.19, 0, 0.019];

[t, y] = ode45(@system, [t0, t1], initials);

total_population = y(:, 1) + y(:, 2);

moi = y(:, 3) ./ y(:, 1);

figure(1);
plot(t, y, t, total_population, t, moi);
ylim([0, 0.3]);
legend(["E: sensitive bacteria", "I: infected bacteria", "C: phage concentration", "Total bacteria", "MOI"]);

function dydt = system(t_system, y)
    alpha_gomez = 1.63;
    a_gomez = 0.4428;
    beta_gomez = 3.81;
    v_gomez = 1.002;
    gamma_gomez = 1.032 * 10^(-2);
    
    %%copyed

    T = @(t) 18; %temperature
    a = 0.11; %scale factor
    lambda = 0.08; %bacterial inactivation rate
    gamma = 0.001; %bacterial activation rate
    
    %mu_g parameters
    alph = 1.265 * 10^(-4); %scale factor
    T_0 = 2.85; %minimum temperature for growth
    K = inf; %carrying capacity. K = inf for exponential growth.
    
    L = @(t) 400; %light intensity
    b = 100; %light intensity where the effect on decay is at 50%
    
    relative_humidity = @(t) 78; %going by the dew point calculation, this isn't supposed to be a percentage
    %dew point calculation from supplemental material
    A_1 = 17.625;
    B_1 = 243.04;
    dew_point = @(T, rh) B_1 * (log(rh / 100) + A_1 * T / (B_1 + T)) / (A_1 - log(rh / 100) - A_1 * T / (B_1 + T)); 
    %dew point depression
    d = @(t) T(t) - dew_point(T(t), relative_humidity(t));
    %dew point depression threshold
    T_DPD = 2;
    
    %I think these need to be computed from  F_lambda & F_gamma?
    q_t_g = 1/2;
    c_t_d = 1/2;
    t_g = 5;
    t_d = 5;
    
    %% the actual model
    
    %bacterial growth rate
    mu_g = @(T, y) alph * (T - T_0)^2 * (1 - gamma / K);
    
    F_gamma = @(t) 1 - (1 - q_t_g) * exp(-gamma * (t - t_g));
    
    %makes heaviside act like the correct type of comparison
    oldparam = sympref('HeavisideAtOrigin', 1);

    mu = @(t, y) mu_g(T(t), y) * F_gamma(t) * heaviside(-d(t) + T_DPD);
    
    F_lambda = @(t) c_t_d * exp(-lambda * (t - t_d));
    
    k = @(t) (1 + L(t) / (b + L(t))) * d(t) * a * F_lambda(t) * heaviside(d(t) - T_DPD);

    %restores your original heaviside origin setting
    sympref('HeavisideAtOrigin', oldparam);

    %%

    dydt = zeros(3, 1);

    %If infected bacteria can't reproduce:
    dydt(1) = (mu(t_system, y(1)) - k(t_system)) * y(1) - beta_gomez * y(1) * y(3) / (a_gomez + y(3));
    %If infected bacteria can reproduce:
    %dydt(1) = (mu(t_system, y(1)) - k(t_system)) * (y(1) + y(2)) - beta_gomez * y(1) * y(3) / (a_gomez + y(3));

    dydt(2) = beta_gomez * y(1) * y(3) / (a_gomez + y(3)) - v_gomez * y(2);
    dydt(3) = alpha_gomez * y(2) - gamma_gomez * y(3);
end