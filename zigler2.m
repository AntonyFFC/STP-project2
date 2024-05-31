K0 = 6.4;
T0 = 5.0;
T1 = 2.07;
T2 = 4.6;
Tp = 0.5;
% Definicja transmitancji obiektu
s = tf('s');
G_s = (K0 * exp(-T0 * s)) / ((T1 * s + 1) * (T2 * s + 1));

% Czas symulacji
t = 0:0.1:1500;

% Eksperymentalne dobieranie wartości krytycznej wzmocnienia
K_critical = 0.33136;

G_with_K = feedback(G_s * K_critical, 1);
[y, t] = step(G_with_K, t);
[pks, locs] = findpeaks(y, t);
T_k = mean(diff(locs));
fprintf('Okres krytyczny T_k = %.4f\n', T_k);
    
% Wykres odpowiedzi skokowej
figure;
plot(t, y, 'b');
title('Odpowiedź skokowa');
xlabel('Czas (s)');
ylabel('y(t)');
print ('pics/odp_skok_zigler.png', '-dpng', '-r400')
grid on;

% Ustawienie układu osi na wykresie sygnału wejściowego
axis([0 1500 -0.5 2]);
