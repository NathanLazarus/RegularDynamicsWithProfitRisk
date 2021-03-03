function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = stage1model(u)

syms DELTA ALFA BETTA G LAMBDAZ ETA MU LAMBDAP LAMBDAphi LAMBDAphi_lagged phi_reg_intercept LAMBDAphi_theta LAMBDAphi_lagtheta GAMA SIGM
syms c cp l lp k kp Z Zp Z1 Z1p Z2 Z2p Z3 Z3p P_lag P_lagp P Pp phi_tminus2 phi_tminus2p phi_tminus1 phi_tminus1p phihat phihatp % riskless_r_ riskless_rp_ risky_rp_ risky_r_

if ~exist('u','var')
    u = (1/(1 - SIGM)) * c^(1-SIGM) * (1 - (GAMA/(1+ETA)) * l^(1+ETA));
end

up = subs(u,[c l],[cp lp]);
dupdcp = jacobian(up,cp);
dudc = jacobian(u,c);
%dudl = jacobian(u,l);

f1 = c + G*kp - (1-DELTA) * k - y_func(k, l, Z, ALFA, phi_tminus1p);
f2 = dudc - BETTA * subs(dupdcp, cp, G*cp) * big_R(kp, lp, Pp, Zp, ALFA, DELTA, phihatp);
f3 = laborsupply(u) + w_func(k, l, P, Z, ALFA, phihat);
f4 = Zp - Z^LAMBDAZ;
f5 = Pp - P_func(G, P, LAMBDAP, Zp, Z, Z1, Z2, Z3, MU, 0) - log(1 + exp(-100*  (P_func(G, P, LAMBDAP, Zp, Z, Z1, Z2, Z3, MU, 0) - 1)))/100;
f6 = Z1p - Z;
f7 = Z2p - Z1;
f8 = Z3p - Z2;
% f9 = riskless_r_ - dudc/(BETTA * subs(dupdcp,cp,G*cp));
% f10 = dudc - BETTA * subs(dupdcp, cp, G*cp) * G * (stockp + ((Pp-1)/Pp)*y_func(kp, lp, Zp, ALFA, phi two shocks in the future))/stock;
% f11 = risky_r_ - big_R(k,l,P,Z,ALFA,DELTA);
f_phi_AR = phi_tminus1p - phi_transition(phi_tminus1, phi_tminus2, P, P_lag, LAMBDAphi, LAMBDAphi_lagged, phi_reg_intercept, LAMBDAphi_theta, LAMBDAphi_lagtheta, 0);
f_phihat = phihat - phi_transition(phi_tminus1, phi_tminus2, P, P_lag, LAMBDAphi, LAMBDAphi_lagged, phi_reg_intercept, LAMBDAphi_theta, LAMBDAphi_lagtheta, 0);
f_phi_lag = phi_tminus2p - phi_tminus1;
f_P_lag = P_lagp - P;

f = [f1;f2;f3;f4;f5;f6;f7;f8;f_phi_AR;f_phihat;f_phi_lag;f_P_lag];

x = [k Z Z1 Z2 Z3 P phi_tminus1 phi_tminus2 P_lag];
y = [l c phihat];
xp = [kp Zp Z1p Z2p Z3p Pp phi_tminus1p phi_tminus2p P_lagp];
yp = [lp cp phihatp];

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);