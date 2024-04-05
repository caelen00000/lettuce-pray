clearvars;
close all;

%% parameters

%parameters from Table 1
%stripped-down model for testing. Only decay conditions needed.
t0 = 0;
t1 = 150;
y0 = 1000000;

T = @(t) 25; %temperature
a = 0.019; %scale factor
lambda = 0.08; %bacterial inactivation rate

L = @(t) 400; %light intensity
b = 70; %light intensity where the effect on decay is at 50%

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
c_t_d = 1/2;
t_d = 0;

%% the actual model

F_lambda = @(t) c_t_d * exp(-lambda * (t - t_d));

k = @(t) (1 + L(t) / (b + L(t))) * d(t) * a * F_lambda(t);

dydt = @(t, y) -k(t) * y;

%% simulation

[t, y] = ode45(dydt, [t0, t1], y0);
figure(1);
semilogy(t, y);
xlabel("time (H)");
ylabel("Log10 population");

%% debug

%{
interval = 1;
ts = 0:interval:t1;
ks = zeros((t1-t0)/interval+1, 1);
for i = t0:interval:t1
    ks(i + 1) = k(i + 1);
end

figure(2);
%plot(t, F_gamma(t), t, F_lambda(t));
plot(ts, ks);
legend("k")
%plot(t, d(t));
%}