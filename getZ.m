function Z=getZ(kx,zmin,zmax,Xm,d,km,k)
N=100000;
dz=(zmax-zmin)/N;
z=zmin:dz:zmax;
F2Z=k*Xm*km*cos(km*(z-d));
Z=interp1(F2Z,z,kx);