function E = getField(X,Z,A,kx,k)
    dkx=kx(2)-kx(1);
    E=zeros(size(Z));
    for i=1:size(kx,2)
        kz=(k^2-kx(i)^2)^0.5;
        E=E+A(i)*exp(1i*(kx(i)*X+kz*Z))*dkx;
    end
end