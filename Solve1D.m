%设定，国际单位制
Nlen=1;
Zmax=0.03;
Zmin=0.01;
Xmin=-600e-6;
Xmax=600e-6;
N=800;
Xm=100e-6;
lambda=632.8e-9;
lambdaM=0.007;
d=10e-3;
dx=(Xmax-Xmin)/N;
dz=(Zmax-Zmin)/N;
z=Zmin:dz:Zmax;
x=Xmin:dx:Xmax;
[X,Z]=meshgrid(x,z);
k=2*pi/lambda;
km=2*pi/lambdaM;
%dkx=2*k/N;
%kx=-k+dkx*0.5:dkx:k;
kx=getKx(Xm,km,k,N);
za=d+lambdaM*(0:0.5:Nlen);
kxf=(kx(1:end-1)+kx(2:end))/2;
Fii=zeros(size(za,2)-1,size(kx,2));
PPi=Fii;
FiPi=PPi;
dkx=kx(2)-kx(1);
for i=1:size(za,2)-1
    Zx=getZ(kx,za(i),za(i+1),Xm,d,km,k);
    Zxf=getZ(kxf,za(i),za(i+1),Xm,d,km,k);
    PP=(1-isnan(Zx))./(k*Xm*km*km*sin(km*(Zx-d)));
    PP(isnan(PP))=0;
    PPi(i,:)=(abs(PP)).^(0.5);
    Fi=-Xm*sin(km*(Zxf-d))+Xm*km*Zxf.*cos(km*(Zxf-d));
    FiP=-Xm*sin(km*(Zx-d))+Xm*km*Zx.*cos(km*(Zx-d));
    Fi(isnan(Fi))=0;
    FiP(isnan(FiP))=0;
    
    FiPi(i,:)=FiP;
    for j=1:size(kxf,2)
        Fii(i,j+1)=Fii(i,j)+Fi(j)*dkx;
    end
end
%Fii(2,:)=Fii(2,:)-Fii(2,1)+Fii(1,1);
PPi(isnan(PPi))=0;
Fii(isnan(Fii))=0;
Pi=PPi.^0.5;
A=PPi.*exp(1i*Fii);
AM=(exp(1i*3*pi*(0:0.5:Nlen-0.1))*A);

AMM=AM.*(abs(kx)<max(abs(kx))*0.9);
E=getField(X,Z,AMM,kx,k);
mesh(Z,X,abs(E)-max(abs(E(:)))/2);
view(2);
hold on;

%scatter3(z,Xm*sin((z-d)*km),ones(size(z))*max(abs(E(:))),'r.')