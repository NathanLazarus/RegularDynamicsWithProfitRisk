% get_coefs_SchmittGrohe_Uribe

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR; Z1 = ZSTAR; Z2 = ZSTAR; Z3 = ZSTAR; P = PSTAR;
kp=k; cp=c; lp=l; Zp=Z; Z1p = Z; Z2p = Z; Z3p = Z; Pp = P; riskless_r_ = RSTAR;
stockSTAR = G*((PSTAR-1)/PSTAR)*y_func(KSTAR,LSTAR,ZSTAR,ALFA,PHISTAR)/((1+RSTAR) - G); stock = stockSTAR; stockp = stockSTAR;
risky_r_ = RSTAR; risky_rp_ = RSTAR; lhat = LSTAR; lhatp = lhat;

Psi1_tminus1 = PSISTARS(1);  Psi2_tminus1 = PSISTARS(2);  Psi3_tminus1 = PSISTARS(3);  Psi4_tminus1 = PSISTARS(4);  Psi5_tminus1 = PSISTARS(5);
Psi1_tminus1p = PSISTARS(1); Psi2_tminus1p = PSISTARS(2); Psi3_tminus1p = PSISTARS(3); Psi4_tminus1p = PSISTARS(4); Psi5_tminus1p = PSISTARS(5);

Psi1 = PSISTARS(1);  Psi2 = PSISTARS(2);  Psi3 = PSISTARS(3);  Psi4 = PSISTARS(4);  Psi5 = PSISTARS(5);
Psi1p = PSISTARS(1); Psi2p = PSISTARS(2); Psi3p = PSISTARS(3); Psi4p = PSISTARS(4); Psi5p = PSISTARS(5);

eta_without_Psis = [0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 sigma_P/sigma_Z 0]; %Matrix defining driving force
eta = zeros([6 + n_Psis 3 + n_Psis]);
eta(1:6,1:3) = eta_without_Psis;
eta(7:(6 + n_Psis), 4:(3 + n_Psis)) = Psi_vcov/sigma_Z;

flatten = @(A) A(:);

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = stage1model(u);

approx = order;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

[nstate,~] = size(hx);

if order >= 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);

    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);

    dec_l_hat = [LSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
%     dec_phi_hat = [PHISTAR,gx(3,:),1/2*flatten(gxx(3,:,:))',1/2*gss(3)];

else
    dec_l_hat = [LSTAR,gx(1,:),zeros([1 nstate^2+1])];
%     dec_phi_hat = [PHISTAR,gx(3,:),zeros([1 nstate^2+1])];

end


eta = [eta;zeros([1 size(eta,2)])]; %adding lhat state variable but removing P_lag


[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = stage2model(u, dec_l_hat, decision_func_to_use, ssvals_stage1);

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

[nstate,~] = size(hx);

if order >= 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx);

    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);

    dec_k=[KSTAR,hx(1,:),1/2*flatten(hxx(1,:,:))',1/2*hss(1)];
    dec_c=[CSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
    dec_stock = [stockSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];

else
    dec_k=[KSTAR,hx(1,:),zeros([1 nstate^2+1])]; 
    dec_c=[CSTAR,gx(1,:),zeros([1 nstate^2+1])];
    dec_stock = [stockSTAR,gx(2,:),zeros([1 nstate^2+1])];

end