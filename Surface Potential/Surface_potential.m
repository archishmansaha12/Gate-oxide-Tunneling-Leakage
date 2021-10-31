clc;
clear all

%% --------parameters--------

nm = 1e-9;  %define 1 nm length
q = 1.6e-19; %unit charge (1e)
m0 = 9.1e-31;

%Vg = 0.03;  %Gate voltage (INPUT)

Ea = 4.1*q; %Electron affinity of Si substrate = 4.1eV
Eg = 1.1*q; %Bandgap of Si = 1.1eV
offset = 3.25*q; %Ec offset between Si and SiO2 = 3.25eV
Eg_ox = 9*q;  %Bandgap of SiO2 = 9eV

t_oxide = [1.5*nm 3*nm];  %oxide thickness (INPUT)
N_sub = 5e18*1e6; %enter -ve for NMOS and +ve for PMOS (INPUT)
rho_fix = 0;%2*1e12*1e4; %fixed oxide charge density (INPUT)
D_it = 0;%3*1e12*1e4; %interface trap density (INPUT)

eps_0 = 8.85e-12; %permittivity of free space
eps_sub = 11.6*eps_0; %permittivity of Si substrate
eps_ox = 3.9*eps_0; %permittivity of SiO2

ni = 1e10*1e6; %intrinsic carrier concentration
p_sub = ni;
n_sub = ni;

kbT = 25.6e-3*q; 

Nc = ni*exp(Eg/(2*kbT)); %Nc of Si

p_sub = p_sub + N_sub; %p = ni + Na
phi_b = (kbT/q)*log(p_sub/ni); %calculates phi_b
n_sub = ni^2/p_sub;   %n = ni^2/p
phi_g = Ea;   %WF of gate (band edge)    

Ecf =  kbT*log(Nc/n_sub); %Gap between CB and Fermi level
phi_s = (Ea+Ecf);  %WF of substrate



Vfb = (phi_g-phi_s)/q; %Flat band voltage

Vg = 0.01:0.01:3;
no_points = length(Vg);
psi_poisson = zeros(length(t_oxide),no_points);
psi_schro_poiss = 0.01:0.001:1.8;
no_points2 = length(psi_schro_poiss);
Vg_schro_poiss = zeros(length(t_oxide),no_points2);

Section1 = 1

for k=1:length(t_oxide)

t_ox = t_oxide(k);    

%% --------calculation of surface potential-----------poisson--------

for j = 1:no_points

    Vg1 = Vg(j)
    [psi_s, F_ox] = poisson(m0,q,eps_sub,eps_ox,t_ox,N_sub,Eg,ni,Vg1,phi_g); 
    psi_poisson(k,j) = psi_s;

end

%% -------------calculation of surface potential---------schrodinger-poisson

for i=1:no_points2

    psi1 = psi_schro_poiss(i)
    [Vg2, F_ox] = schrodinger_poisson(m0,q,eps_sub,eps_ox,t_ox,N_sub,Eg,ni,psi1,phi_g);
    Vg_schro_poiss(k,i) = Vg2;
    
end

end

%% -------------------plotting---------------------------

figure;
plot(Vg,psi_poisson(1,:),'k:','linewidth',1.8);
hold on;
plot(Vg_schro_poiss(1,:),psi_schro_poiss,'k','linewidth',1.1);
hold on;
plot(Vg,psi_poisson(2,:),'k--','linewidth',2);
hold on;
plot(Vg_schro_poiss(2,:),psi_schro_poiss,'k','linewidth',2);
xlim([0 3]);
grid on;
xlabel('Gate Voltage   V_{G} (in V)');
ylabel('Surface Potential   \phi_s (V)');
legend('t_{ox}=1.5nm Poisson solver','t_{ox}=1.5nm Schrodinger-Poisson solver','t_{ox}=3.0nm Poisson solver','t_{ox}=3.0nm Schrodinger-Poisson solver','Location','best');
title('Surface Potential vs Gate Voltage');
