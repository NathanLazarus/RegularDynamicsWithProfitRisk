% model1.M
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = stage2model(u, dec_l_hat, decision_func_to_use, ssvals)

syms DELTA ALFA BETTA G LAMBDAZ ETA MU LAMBDAP LAMBDAphi GAMA SIGM sigma_Z sigma_P
syms c cp k kp Z Zp Z1 Z1p Z2 Z2p Z3 Z3p P Pp stock stockp phi phip lhat lhatp l % riskless_r_ riskless_rp_ risky_rp_ risky_r_

if ~exist('u','var')
    u = (1/(1 - SIGM)) * c^(1-SIGM) * (1 - (GAMA/(1+ETA)) * l^(1+ETA));
end

% lhatp = f(kp, Zp, Z, Z1, Z2, Z3, Pp, phi);
% rhatp = g(kp, Zp, Z, Z1, Z2, Z3, Pp, phi);
% rhatp = phi^LAMBDAphi .* (1 ./ Pp) .* ALFA .* Zp .* kp^(ALFA - 1) * lhatp ^ (1 - ALFA);
rhatp = big_R(kp, lhatp, Pp, Zp, ALFA, DELTA, phihatp);

up = subs(u,[c l],[cp lhatp]);
dupdcp = jacobian(up,cp);
u_t = subs(u,l,lhat);
dudc = jacobian(u_t,c);
%dudl = jacobian(u,l);

f1 = c + G*kp - (1-DELTA) * k - y_func(k, lhat, Z, ALFA, phi);
f2 = dudc - BETTA * subs(dupdcp, cp, G*cp) * rhatp;
% f3 = laborsupply(u) + w_func(k, l, P, Z, ALFA, phihat);
f4 = Zp - Z^LAMBDAZ;
f5 = Pp - P_func(G, P, LAMBDAP, Zp, Z, Z1, Z2, Z3, MU, 0) - log(1 + exp(-100*  (P_func(G, P, LAMBDAP, Zp, Z, Z1, Z2, Z3, MU, 0) - 1)))/100;
f6 = Z1p - Z;
f7 = Z2p - Z1;
f8 = Z3p - Z2;
% f9 = riskless_r_ - dudc/(BETTA * subs(dupdcp,cp,G*cp));
f10 = dudc - BETTA * subs(dupdcp, cp, G*cp) * G * (stockp + ((Pp-1)/Pp)*y_func(kp, lhatp, Zp, ALFA, phip))/stock;
% f11 = risky_r_ - big_R(k,l,P,Z,ALFA,DELTA);
f_phi_AR = phip - phi_transition(phi, Pp);
f_phihat = phihatp - phi_transition(phi, Pp);
f_lhat = lhatp - decision_func_to_use(dec_l_hat, [kp, Zp, Z, Z1, Z2, Pp, phihatp], ssvals, sigma_Z, sigma_P);

f = [f1;f2;f4;f5;f6;f7;f8;f10;f_phi_AR;f_lhat];

x = [k Z Z1 Z2 Z3 P phi lhat];
y = [c stock phihat];
xp = [kp Zp Z1p Z2p Z3p Pp phip lhatp];
yp = [cp stockp phihatp];

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);