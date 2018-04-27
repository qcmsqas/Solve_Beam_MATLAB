%设定，国际单位制
Nlen=3;
d=8e-3;
Zmax=d+2e-3;
Zmin=d-2e-3;
Xmin=-100e-6;
Xmax=100e-6;
N=400;
Xm=100e-6;
lambda=632.8e-9;
lambdaM=0.07;
dx=(Xmax-Xmin)/N;
dz=(Zmax-Zmin)/N;
z=Zmin:dz:Zmax;
x=Xmin:dx:Xmax;
Lu=1.5e-3;
a=Xm/Lu/Lu;
k=2*pi/lambda;
km=2*pi/lambdaM;
dzz=Lu/N;
zz=d-Lu:dzz:d+Lu;
kx=(zz-d)*2*k*a;
[X,Z]=meshgrid(x,z);
%dkx=2*k/N;
%kx=-k+dkx*0.5:dkx:k;
%kx=getKx(Xm,km,k,N);
%za=d+lambdaM*(0:0.25:Nlen);
kxf=(kx(1:end-1)+kx(2:end));
%Fii=zeros(size(za,2)-1,size(kx,2));
%PPi=Fii;
%FiPi=PPi;
%for i=1:size(za,2)-1
%    Zx=getZ(kx,za(i),za(i+1),Xm,d,km,k);
%    Zxf=getZ(kxf,za(i),za(i+1),Xm,d,km,k);
%    PPi(i,:)=1./(k*Xm*km*km*sin(km*(Zx-d)));
%    Fi=Xm*sin(km*(Zxf-d))-Xm*km*Zxf.*cos(km*(Zxf-d));
%    FiP=Xm*sin(km*(Zx-d))-Xm*km*Zx.*cos(km*(Zx-d));
%    FiPi(i,:)=FiP;
%    for j=1:size(kxf,2)
%        Fii(i,j+1)=Fii(i,j)+Fi(j);
%    end
%end
%PPi(isnan(PPi))=0;
%Fii(isnan(Fii))=0;
%Pi=PPi.^0.5;
Fi=-kxf.*kxf/(4*k*k*a);
Fii=kx.^3/(12*k*k*a)+kx.^2*d/k/2;
%dkx=kx(2)-kx(1);
%for j=1:size(kxf,2)
%    Fii(j+1)=Fii(j)+Fi(j)*dkx;
%end
P=ones(size(kx));
E = getField(X,Z,P.*exp(1i*Fii),kx,k);
mesh(Z,X,abs(E));
view(2)