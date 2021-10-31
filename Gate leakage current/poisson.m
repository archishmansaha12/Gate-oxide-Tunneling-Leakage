function [psi_s, F_ox] = poisson(m0,q,eps_sub,eps_ox,t_ox,N_dop,Eg,ni,Vg,phi_g)

%% ---------------Parameters-----------------------------

m_l = 0.35*m0;
n_vj = 1;
m_dj = 0.025*m0;
m_yj = 0.98*m0;
h = 6.626e-34;
hbar = h/2/pi;
kbT = 25.6e-3*q;

E_cf = kbT*log(N_dop/ni);
phi_s = phi_g+Eg/2+E_cf;
Vfb = (phi_g-phi_s)/q;
p_sub = ni+N_dop;
n_sub = ni^2/p_sub;

%% --------- solve Poisson equation-----------------------------

N1_iter1 = 1000; %No of iterations
psi_s = sign(Vg-Vfb)*0.1; %initial value of surface potential
C_ox = eps_ox/t_ox;

for i1 = 1:N1_iter1  %Newton Rhapson method to solve Poisson equation
    
    Qs_LF = -sign(Vg-Vfb)*sqrt(2*eps_sub*kbT)*sqrt(abs(p_sub*(exp(-q*psi_s/kbT) + q*psi_s/kbT -1)+ n_sub*(exp(q*psi_s/kbT) - q*psi_s/kbT -1)));
    a = Vfb - Vg - Qs_LF/C_ox + psi_s;
    del_a = 1 + sign(Vg-Vfb)*(1/C_ox)*sqrt(eps_sub*kbT/2)*(1/sqrt(abs(p_sub*(exp(-q*psi_s/kbT) + q*psi_s/kbT -1)+ n_sub*(exp(q*psi_s/kbT) - q*psi_s/kbT -1))))*(p_sub*q/kbT*(1-exp(-q*psi_s/kbT)) + n_sub*q/kbT*(exp(q*psi_s/kbT)-1));
    del_psi_s = -a/del_a;
    psi_s = psi_s + del_psi_s;
    
end

psi_s_poisson = psi_s;

% %% -------------initializing for Schrodinger-Poisson equation------------
% 
% E_fn = -q*(Eg/(2*q) - psi_s + E_cf/q);
% N_dep = sqrt(2*eps_sub*psi_s*N_dop/q);
% % if psi_s > 2*E_cf/q
% % N_inv = C_ox*(psi_s-2*E_cf/q);%N_dep;
% % else
% N_inv = N_dep;    
% % end
% F_dep = q*N_dep/eps_sub;
% F_inv = q*N_inv/eps_sub;
% F_s = F_dep+F_inv;
% 
% i_max=3; j_max=2;
% 
% N1_iter2 = 5e3;
% 
% for i2=1:N1_iter2
%    
%     i = (1:1:i_max)'.*ones(1,j_max);
%     
%     E_ij_dep = (hbar^2/2/m_yj)^(1/3)*(3/2*pi*q*F_s*(i-0.25)).^(2/3);
%     b = ((12*m_l*q^2/eps_sub/h^2)*(N_dep+11/32*N_inv))^(1/3);
%     Z0 = 3/b;
%     
%     %%% calculating energies of levels
%     E_ij = E_ij_dep - q^2*F_dep*F_inv*Z0^2/4./E_ij_dep - 4*E_ij_dep.^2/15/q/F_dep/Z0 + q*F_inv*Z0;
%     E_11 = (1.5)^(5/3)*(q^2*hbar/sqrt(m_l)/eps_sub)^(2/3)*(N_dep+55/96*N_inv)*(N_dep+11/32*N_inv)^(-1/3);
%     E_ij(1,1) = E_11;
%     
%     N_ij = n_vj*m_dj*kbT/pi/hbar^2*log(1+exp((E_fn-E_ij)/kbT));
%     
%     Z_ij = 2/3*E_ij/q/F_s;
%     
%     N_inv = sum(sum(N_ij)); %%% inversion charge
%     
%     Z_av = sum(sum(N_ij.*Z_ij))/N_inv;
%     psi_dep = psi_s - kbT/q - q*N_inv*Z_av/eps_sub;
%     N_dep = sqrt(2*eps_sub*psi_dep*N_dop/q);
%     
%     F_dep = q*N_dep/eps_sub;
%     F_inv = q*N_inv/eps_sub;
%     F_s = F_dep+F_inv;
%     
%     
% %     if Vg - Vfb - V_ox > 0
% %     psi_s = Vg - Vfb - V_ox;
% %     end
% %     E_fn = -q*(Eg/(2*q) - psi_s + E_cf/q);
%     
% end    

V_ox = Vg - Vfb - psi_s;
F_ox = V_ox/t_ox;
% F_ox = F_s*eps_sub/eps_ox;
% V_ox = F_ox*t_ox;

end