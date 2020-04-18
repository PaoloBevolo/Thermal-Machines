%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Model of Thermal Machine MACI - course @ UPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PROPIEDADES

%Fluido motor aire
R  = 8.31 %[J/mol K]
disp('[J/mol k]')
Cp = 28.14 %[J/mol K]
Cv = Cp - R %[J/mol K]
M  = 28.9 %[grs/mol]
M  = M/1000; %[kg/mol]
k  = Cp/Cv %[]
ENT = 0.179 %[J/ mol·K]

%Geometria cilindro
Vcu = 400 %[cm3]
Vcu = Vcu / 1000000; %[m3]
Vcc = 40 %[cm3]
Vcc = Vcc / 1000000; %[m3]
rel = (Vcc + Vcu) / Vcc %[]
D = 10 %[cm]
D = D / 100; %[m]
CA = Vcu / (pi*D^2/4) %[m2]
S = pi*D^2/4 %[m2]
RPM = 5000 %[RPM]
RPM = RPM/60*2*pi; %[rad/s]
l = 1.8 * CA %[m]
r = CA / 2 %[m]
lambda = r / l %[]

%Propiedades combustible
dos = 14.7 %[]
PCI = 42 %[MJ/kg]
PCI = PCI * 1000000; %[J/kg]

%Variables dinamica
alphamin = 0 %[rad]
alphamax = 2*pi %[rad]
N = 1000 %[]
alpha = linspace(alphamin, alphamax, N); %[rad]
alpha_grad = alpha*180/pi;

%Distribucion
ini_comp = 0 %[º]
ini_comp = ini_comp/180*pi; %[rad]
ini_exp = 180 %[º]
ini_exp = ini_exp/180*pi; %[rad]
ini_comb = 180 %[º]
ini_comb = ini_comb/180*pi; %[rad]

ind_comp = find(alpha > ini_comp & alpha < ini_comb);
ind_exp = find(alpha >= ini_comb);
ind_comb = ind_exp(1);
%ind_exp = ind_exp(2:end);

%Modelo de combustion
IC = 35 %[º]
IC = IC/180*pi; %[rad]
AC = 90 %[º]
AC = AC/180*pi; %[rad]
a_wiebe = 7
m_wiebe = 3
alpha_0 = ini_comb - IC;
alpha_fin = alpha_0 + AC;


%% Cinematica
x = r * (1 - cos(alpha)) + l * (1 - sqrt(1 - lambda^2*sin(alpha).^2));
v = r * (sin(alpha) + lambda * sin(alpha).*cos(alpha) / sqrt(1 - lambda^2*sin(alpha).^2)) * RPM;
a = RPM^2 * r * (cos(alpha) + lambda * cos(2*alpha));
V = Vcc + (CA - x) * S; %[m^3]

%Plot cinematica
figure
plot(alpha_grad,x)
xlabel('\alpha [º]')
ylabel('Distancia [m]')
grid on
set(gcf,'Color',[1 1 1]);

figure
plot(alpha_grad,v)
xlabel('\alpha [º]')
ylabel('Velocidad [m/s]')
grid on
set(gcf,'Color',[1 1 1]);

figure
plot(alpha_grad,a)
xlabel('\alpha [º]')
ylabel('Aceleracion [m/s^2]')
grid on
set(gcf,'Color',[1 1 1]);

figure
plot(alpha_grad,V*1000000)
xlabel('\alpha [º]')
ylabel('Volumen [cm^3]')
grid on
set(gcf,'Color',[1 1 1]);


%% CONDICIONES INICIALES
T1 = 15 %[º]
T1 = T1 + 273; %[K]
P1 = 1 %[atm]
P1 = P1 * 100000; %[Pa]
V1 = V(1); %[m^3]
T4 = T1;
P4 = P1;
V4 = V1;

%% Termodinamica
n = P1 * (Vcc + Vcu) / R / T1 %[mol]
S1 = ENT * n
S4 = S1
masa_air = n * M * 1000 %[grs]
masa_fuel = masa_air / dos %[grs]
masa_fuel = masa_fuel / 1000; %[kg]
delta_Q = masa_fuel * PCI %[J]
delta_T = delta_Q / Cv / n %[K]

%No combustion
P_nc = P1 * (V(1) ./ V).^k; %[Pa]
T_nc = P_nc .* V / n / R; %[K]
S_nc(1:N) = S1; %[J/T]

%Combustion ideal
P_teor = P_nc; %[Pa]
T_teor = T_nc; %[K]
S_teor = S_nc; %[J/T]
V2 = V(ind_comb)
V3 = V2
T2 = T_teor(ind_comb)
T3 = T2 + delta_T;
P2 = P_teor(ind_comb);
P3 = n * R * T3 / V(ind_comb)
S2 = S1; %[J/K]
S3 = S2 + Cv * log(T3/T2) + R * log(V3/V2); %[J/K]
P_teor(ind_exp) = P3 * (V(ind_comb) ./ V(ind_exp)).^k; %[Pa]
T_teor(ind_exp) = P_teor(ind_exp) .* V(ind_exp) / n / R; %[K]
S_teor(ind_exp) = S3; %[J/K]
P_teor(end) = P4;
T_teor(end) = T4;
S_teor(end) = S4;

P_isot = P_teor;
P_isot(ind_exp) = n * R * T3 ./ V(ind_exp);


%% Modelo de combustion
alpha_bar = (alpha - alpha_0) / AC;
alpha_ind = find(alpha_bar >=0 & alpha_bar <= 1);
alpha_bar = alpha_bar(alpha_ind);
y_wiebe = 1 - 1 ./ exp(a_wiebe * alpha_bar.^(m_wiebe + 1));
delta_y_wiebe = [diff(y_wiebe) 0];
delta_Q_wiebe = zeros(N,1);
delta_Q_wiebe(alpha_ind) = delta_y_wiebe * delta_Q;

%Plot modelo de combustion
figure
plot(alpha_bar,y_wiebe)
line(alpha_bar,N/(alphamax - alphamin)*delta_y_wiebe,'Color','r')
xlabel('\alpha_{bar} []')
ylabel('Y [%]')
grid on
set(gcf,'Color',[1 1 1]);

%Wiebe
P_real = P_nc; %[Pa]
T_real = T_nc; %[K]
S_real = S_nc; %[J/T]

%Nusselt
P_real_HT = P_nc; %[Pa]
T_real_HT = T_nc; %[K]
S_real_HT = S_nc; %[J/T]


delta_V = [0 diff(V)];

%Loop principal
for i=2:N,

    P_real(i) = P_real(i-1) + (delta_Q_wiebe(i) - Cp / R * P_real(i-1) * delta_V(i)) * R / Cv / V(i);
    T_real(i) = P_real(i) * V(i) / n / R;
    %S_real(i) = S_real(i-1) + delta_Q_wiebe(i) / T_real(i-1);
    S_real(i) = S_real(i-1) + Cv * log(T_real(i)/T_real(i-1)) + R * log(V(i)/V(i-1));

    %Extracción de calor
    perdidas = 1500 * pi * D * CA / 2 * (T_real_HT(i-1) - 90 + 273) * (2 * pi) / RPM / N ;
    P_real_HT(i) = P_real_HT(i-1) + (delta_Q_wiebe(i) - perdidas - Cp / R * P_real_HT(i-1) * delta_V(i)) * R / Cv / V(i);
    T_real_HT(i) = P_real_HT(i) * V(i) / n / R;
    S_real_HT(i) = S_real_HT(i-1) + Cv * log(T_real_HT(i)/T_real_HT(i-1)) + R * log(V(i)/V(i-1));

end


%Plot termodinamica
figure
h_presion = plot(alpha_grad,P_nc/100000,'LineWidth',2)
line(alpha_grad,P_teor/100000,'Color','r')
line(alpha_grad,P_real/100000,'Color','k')
xlabel('\alpha [º]')
ylabel('Presion [atm]')
grid on
set(gcf,'Color',[1 1 1]);
legend('No combustion', 'Ideal', 'Real');

figure
h_PV = plot(V*1000000,P_nc/100000,'LineWidth',2)
line(V*1000000,P_teor/100000,'Color','r')
line(V*1000000,P_real/100000,'Color','k','LineWidth',2)
line(V(ind_exp)*1000000,P_isot(ind_exp)/100000,'Color','k','LineStyle','-.')
xlabel('Volume [cm^3]')
ylabel('Presion [atm]')
grid on
set(gcf,'Color',[1 1 1]);
legend('No combustion', 'Ideal', 'Real', 'Isot');

figure
h_TS = plot(T_nc - 273,S_nc,'LineWidth',2)
line(T_nc - 273,S_teor,'Color','r')
line(T_nc - 273,S_real,'Color','k')
xlabel('T [º]')
ylabel('S [J/K]')
grid on
set(gcf,'Color',[1 1 1]);
legend('No combustion', 'Ideal', 'Real');



figure
h_presion = plot(alpha_grad,P_real/100000,'LineWidth',2)
line(alpha_grad,P_real_HT/100000,'Color','r')
xlabel('\alpha [º]')
ylabel('Presion [atm]')
grid on
set(gcf,'Color',[1 1 1]);
legend('Wiebe', 'Wiebe + Nusselt');

figure
h_PV = plot(V*1000000,P_real/100000,'LineWidth',2)
line(V*1000000,P_real_HT/100000,'Color','r')
xlabel('Volume [cm^3]')
ylabel('Presion [atm]')
grid on
set(gcf,'Color',[1 1 1]);
legend('Wiebe', 'Wiebe + Nusselt');

figure
h_TS = plot(T_nc - 273,S_real,'LineWidth',2)
line(T_nc - 273,S_real_HT,'Color','r')
xlabel('T [º]')
ylabel('S [J/K]')
grid on
set(gcf,'Color',[1 1 1]);
legend('Wiebe', 'Wiebe + Nusselt');


W_teor = P_teor .* [0 diff(V)]
W_real = P_real .* [0 diff(V)]

nu_teor = sum(W_teor) / delta_Q * 100
nu_real = sum(W_real) / delta_Q * 100
