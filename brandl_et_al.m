clearvars;
close all;

%% parameters

%parameters from 'Model application to EcO157 dynamics..', morning trials
%several of these come from data from real world climate conditions, I just
%used averages or values from earlier in the paper for now. Thats why
%there are so many constant anonymous functions of time.
t0 = 0;
t1 = 150;
y0 = 1000000;

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

%Is this correct? Does this need a different heaviside 0 value?
mu = @(t, y) mu_g(T(t), y) * F_gamma(t) * heaviside(-d(t) + T_DPD);

F_lambda = @(t) c_t_d * exp(-lambda * (t - t_d));

%Is this correct? Does this need a different heaviside 0 value?
k = @(t) (1 + L(t) / (b + L(t))) * d(t) * a * F_lambda(t) * heaviside(d(t) - T_DPD);

dydt = @(t, y) (mu(t, y) - k(t)) * y;

%% simulation

[t, y] = ode45(dydt, [t0, t1], y0);
figure(1);
semilogy(t, y);
xlabel("time (H)");
ylabel("Log10 population");

%% debug

%{
t = 1:100;
ks = zeros(100);
mus = zeros(100);
for i = t 
    ks(i) = k(i);
    mus(i) = mu(i, 0);
end

%figure(2);
%plot(t, F_gamma(t), t, F_lambda(t));
%plot(t, ks);
%plot(t, mus);
%plot(t, d(t));
%}