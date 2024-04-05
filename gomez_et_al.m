clearvars;
close all;

t0 = 0;
t1 = 4;

initials = [0.19, 0, 0.019];

[t, y] = ode45(@system, [t0, t1], initials);

total_population = y(:, 1) + y(:, 2);

figure(1);
plot(t, y, t, total_population);
legend(["E: sensitive bacteria", "I: infected bacteria", "C: phage concentration", "Total bacteria"]);

function dydt = system(~, y)
    k = 1.11;
    alpha = 1.63;
    a = 0.4428;
    beta = 3.81;
    v = 1.002;
    gamma = 1.032 * 10^(-2);
    N = 3.2;

    dydt = zeros(3, 1);
    dydt(1) = k * y(1) * (1 - (y(1) + y(2)) / N) - beta * y(1) * y(3) / (a + y(3));
    dydt(2) = beta * y(1) * y(3) / (a + y(3)) - v * y(2);
    dydt(3) = alpha * y(2) - gamma * y(3);
end