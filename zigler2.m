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
