% DMC
%clear all;
%zad1
Dmc = struct;
Dmc.czas_symulacji = 100;
Dmc.y_zad = 1;
c = [0.0748, 0.0666];
b = [1, -1.6824, 0.7045];
Dmc.wyu = 0;
Dmc.wyy = 0;
Dmc.lambda = 1;
%horyzonty 
Dmc.D = 69;
Dmc.N = 69;
Dmc.Nu = 69;
Dmc.y_start = 0;
Dmc.u_start = 0;
y(1:12) = Dmc.y_start;
u(1:12) = Dmc.u_start;
s = step(G_z);
s(1) = [];
z(1:40) = 0;
z(41:Dmc.czas_symulacji) = 0.05;
for i=1:Dmc.D-1
    deltaupk(i) = 0;
end


% matrix generation
M = zeros(Dmc.N, Dmc.Nu);
for i=1:Dmc.N
    for j = 1:Dmc.Nu
        if (i>=j)
            M(i,j) = s(i-j+1);
        end
    end
end


% mp matrix
Mp = zeros(Dmc.N, Dmc.D-1);
for i = 1:Dmc.N
    for j = 1:(Dmc.D-1)
        if j + i <= Dmc.D
            Mp(i, j) = s(i+j) - s(j);
        else
            Mp(i,j) = s(Dmc.D) - s(j);
        end
    end
end

regulator = struct();
I = eye(Dmc.Nu);
regulator.K =inv( (M' * M + Dmc.lambda * I)) *M';
regulator.Ku = regulator.K(1,:)*Mp;
regulator.Ke = sum(regulator.K(1,:));

for k = 13:Dmc.czas_symulacji
    %obiekt
    y(k) = c(1)*u(k-11) + c(2)*u(k-12) - b(2)*y(k-1) - b(3)*y(k-2);
    %regulator
    ek = Dmc.y_zad - y(k);
    % prawo regulacji
    deltauk = regulator.Ke*ek - regulator.Ku * deltaupk';
    for n=Dmc.D-1:-1:2
        deltaupk(n) = deltaupk(n-1);
    end

    deltaupk(1) = deltauk;
    u(k) = u(k-1)+deltaupk(1);

    Dmc.wyu(k) = u(k);
    Dmc.wyy(k) = y(k);
end

%wizualizacja

figure(1);
stairs(0:Dmc.czas_symulacji, [Dmc.u_start Dmc.wyu]); 
hold on; 
grid on;
xlabel('czas'); 
ylabel('u');
title('Sygnał sterujący')
figure(2);
stairs(0:Dmc.czas_symulacji, [0 Dmc.y_zad*ones(1, Dmc.czas_symulacji)]);
hold on; 
grid on;
plot(1:Dmc.czas_symulacji, Dmc.wyy);
xlabel('czas');
ylabel('y, y_zad');
title('DMC')

