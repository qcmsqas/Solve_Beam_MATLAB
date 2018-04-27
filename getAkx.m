function getAkx(zmin,zmax,k)
N=1000;
dz=(zmax-zmin)/N;
z=zmin+0.5*dz:dz:zmax;
dk=2*k/N;
kx=-k+dk:dk:k;
[Z,Kx]=meshgrid(z,kx);
