% Definicja parametrów obiektu
K0 = 6.4;
T0 = 5;
T1 = 2.07;
T2 = 4.6;
Tp = 0.5; % Okres próbkowania

% Definicja transmitancji ciągłej
s = tf('s');
G_s = (K0 * exp(-T0 * s)) / ((T1 * s + 1) * (T2 * s + 1));

% Dyskretyzacja transmitancji ciągłej
G_z = c2d(G_s, Tp, 'zoh');

% Krytyczne wzmocnienie
K_k = 0.3218;

% Definicja regulatora P z krytycznym wzmocnieniem
C = pid(K_k);

% Symulacja zamkniętego układu regulacji
T = feedback(C * G_z, 1);

% Odpowiedź skokowa
kk = 200; % Liczba próbek
[y, t] = step(T, kk*Tp);

% Znalezienie okresu oscylacji
[pks, locs] = findpeaks(y, t);

% Obliczenie średniego okresu oscylacji
if length(locs) > 1
    T_k = mean(diff(locs));
    fprintf('Okres krytyczny T_k = %.4f\n', T_k);
else
    warning('Nie udało się wykryć pełnych oscylacji.');
    T_k = NaN;
end

% Wyświetlenie wyników
fprintf('Krytyczne wzmocnienie K_k = %.4f\n', K_k);

% Rysowanie odpowiedzi skokowej
figure;
plot(t, y);
title('Odpowiedź skokowa przy krytycznym wzmocnieniu');
xlabel('Czas (s)');
ylabel('Amplituda');
grid on;