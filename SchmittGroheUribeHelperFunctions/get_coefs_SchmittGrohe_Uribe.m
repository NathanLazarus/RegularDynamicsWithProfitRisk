% get_coefs_SchmittGrohe_Uribe

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR; Z1 = ZSTAR; Z2 = ZSTAR; Z3 = ZSTAR; P = PSTAR;
kp=k; cp=c; lp=l; Zp=Z; Z1p = Z; Z2p = Z; Z3p = Z; Pp = P; riskless_r_ = RSTAR;
stockSTAR = G*((PSTAR-1)/PSTAR)*y_func(KSTAR,LSTAR,ZSTAR,ALFA,PHISTAR)/((1+RSTAR) - G); stock = stockSTAR; stockp = stockSTAR;
risky_r_ = RSTAR; risky_rp_ = RSTAR;

lhat = LSTAR; phi = PHISTAR;  phi_tminus1 = PHISTAR;      phi_tminus2 = PHISTAR;      phihat = PHISTAR; P_lag = PSTAR;

lhatp = lhat; phip = PHISTAR; phi_tminus1p = phi_tminus1; phi_tminus2p = phi_tminus2; phihatp = phihat; P_lagp = PSTAR;


eta     = [0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 sigma_P/sigma_Z 0; 0 0 sigma_phi/sigma_Z; 0 0 0; 0 0 0]; %Matrix defining driving force

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
    dec_phi_hat = [PHISTAR,gx(3,:),1/2*flatten(gxx(3,:,:))',1/2*gss(3)];

else
    dec_l_hat = [LSTAR,gx(1,:),zeros([1 nstate^2+1])];
    dec_phi_hat = [PHISTAR,gx(3,:),zeros([1 nstate^2+1])];

end


eta     = [0 0 0; 1 0 0; 0 0 0; 0 0 0; 0 0 0; 0 sigma_P/sigma_Z 0; 0 0 sigma_phi/sigma_Z; 0 0 0; 0 0 0]; %adding lhat state variable but removing P_lag


[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = stage2model(u, dec_l_hat, dec_phi_hat, decision_func_to_use, ssvals_stage1);

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