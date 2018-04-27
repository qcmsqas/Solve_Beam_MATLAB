function kx=getKx(Xm,km,k,Np)
kxMax=k*Xm*km;
dkx=2*kxMax/Np;
kx=-kxMax+0.5*dkx:dkx:kxMax;