% Punkt 1
K0 = 6.4;
T0 = 5.0;
T1 = 2.07;
T2 = 4.6;
Tp = 0.5;

s = tf('s');
G_s = (K0 * exp(-T0 * s)) / ((T1 * s + 1) * (T2 * s + 1));
disp('Transmitancja ciągła G(s):');
disp(G_s);

G_z = c2d(G_s, Tp, 'zoh');
disp('Transmitancja dyskretna G(z):');
disp(G_z);

% Porównanie odpowiedzi skokowej transmitancji ciągłej i dyskretnej
t_cont = 0:0.1:100;  % Czas dla odpowiedzi ciągłej
t_disc = 0:Tp:100;   % Czas dla odpowiedzi dyskretnej

% Odpowiedź skokowa dla transmitancji ciągłej
[y_cont, t_cont] = step(G_s, t_cont);

% Odpowiedź skokowa dla transmitancji dyskretnej
[y_disc, t_disc] = step(G_z, t_disc);

% Wykresy odpowiedzi skokowej
figure;
plot(t_cont, y_cont, 'b', 'LineWidth', 1.5); hold on;
stairs(t_disc, y_disc, 'r', 'LineWidth', 1.5);
legend('Ciągła odpowiedź skokowa', 'Dyskretna odpowiedź skokowa');
xlabel('Czas (s)');
ylabel('Odpowiedź skokowa');
title('Porównanie odpowiedzi skokowej');
grid on;

% Wzmocnienie statyczne dla transmitancji ciągłej
K_static_cont = dcgain(G_s);
disp('Wzmocnienie statyczne transmitancji ciągłej:');
disp(K_static_cont)

% Wzmocnienie statyczne dla transmitancji dyskretnej
K_static_disc = dcgain(G_z);
disp('Wzmocnienie statyczne transmitancji dyskretnej:');
disp(K_static_disc)

% Punkt 2
% Uzyskanie współczynników transmitancji dyskretnej
[num, den] = tfdata(G_z, 'v');

% Wyświetlenie współczynników licznika i mianownika
disp('Współczynniki licznika:');
disp(num);
disp('Współczynniki mianownika:');
disp(den);

% Wyświetlenie równania różnicowego
disp('Równanie różnicowe:');
fprintf('y(k) = ');

% Współczynniki mianownika (a_i) z wyjątkiem pierwszego elementu
for i = 2:length(den)
    if den(i) ~= 0
        fprintf('%+.5f * y(k-%d) ', -den(i), i-1);
    end
end

% Współczynniki licznika (b_i)
for i = 1:length(num)
    if num(i) ~= 0
        if i == 1
            fprintf('%+.5f * u(k) ', num(i));
        else
            fprintf('%+.5f * u(k-%d) ', num(i), i-1);
        end
    end
end

fprintf('\n');

% Punkt 3
% Parametry krytyczne
K_k = 0.3218;
T_k = 20.4;
disp('Parametry krytyczne:');
fprintf('K_k = %.4f\n', K_k);
fprintf('T_k = %.4f\n', T_k);

% Parametry regulatora PID
K_r = 0.6 * K_k;
T_i = 0.5 * T_k;
T_d = 0.12 * T_k;

disp('Parametry regulatora PID:');
fprintf('K_r = %.4f\n', K_r);
fprintf('T_i = %.4f\n', T_i);
fprintf('T_d = %.4f\n', T_d);


% Obliczenie współczynników dyskretnego regulatora PID
r0 = K_r * (1 + Tp / (2 * T_i) + T_d / Tp);
r1 = K_r * (Tp / (2 * T_i) - 2 * T_d / Tp - 1);
r2 = K_r * T_d / Tp;

disp('Współczynniki regulatora dyskretnego PID:');
fprintf('r0 = %.4f\n', r0);
fprintf('r1 = %.4f\n', r1);
fprintf('r2 = %.4f\n', r2);

kk = 300; % koniec symulacji
u = zeros(1, kk); % sygnał sterujący
y = zeros(1, kk); % sygnał wyjściowy
yzad = zeros(1, kk); % sygnał zadany
yzad(10:kk) = 1; % zadana wartość od 10 próbki
e = zeros(1, kk); % uchyb regulacji

% Współczynniki obiektu
b0 = num(1);
b1 = num(2);
a1 = den(2);
a0 = den(3);

for k = 5:kk
    % Symulacja obiektu
    y(k) = b1 * u(k-3) + b0 * u(k-4) - a1 * y(k-1) - a0 * y(k-2);
    
    % Uchyb regulacji
    e(k) = yzad(k) - y(k);
    
    % Sygnał sterujący regulatora PID
    u(k) = r2 * e(k-2) + r1 * e(k-1) + r0 * e(k) + u(k-1);
end

% Wyniki symulacji
figure; stairs(u);
title('u'); xlabel('k');
figure; stairs(y); hold on; stairs(yzad, ':');
title('yzad, y'); xlabel('k');
ylim([-0.0 1.2]);