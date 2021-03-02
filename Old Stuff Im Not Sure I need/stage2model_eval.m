function [dec_k, dec_c, dec_stock] = stage2model_eval(phi, LHAT, order, fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f)


approx = order;

phip = phi ^ 0.95;

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
    dec_c=[CSTAR,gx(1,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
%     dec_riskless_r=[RSTAR,gx(3,:),1/2*flatten(gxx(3,:,:))',1/2*gss(3)];
    dec_stock = [stockSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
%     dec_risky_r = [RSTAR,gx(4,:),1/2*flatten(gxx(4,:,:))',1/2*gss(4)];

else
    dec_k=[KSTAR,hx(1,:),zeros([1 nstate^2+1])]; 
    dec_c=[CSTAR,gx(1,:),zeros([1 nstate^2+1])];
    dec_stock = [stockSTAR,gx(2,:),zeros([1 nstate^2+1])];

end
