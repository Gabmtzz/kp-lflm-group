function Hsimple(Eg,kx,ky,kz,A,B,C,D)
    h=[Eg+D*(kx^2+ky^2+kz^2) 0 0 0 ;
        0 -A*kx^2-B*(ky^2+kz^2) -C*kx*ky -C*kx*kz ;
        0 -C*kx*ky -A*ky^2-B*(kx^2+kz^2) -C*ky*kz ;
        0 -C*kx*kz -C*ky*kz -A*kz^2-B*(kx^2+ky^2)]
    return h
end
function HKane(kx,ky,kz,Eg,Ep,L,N,M,F)
    h4=[Eg+(1+2*F)*(kx^2+ky^2+kz^2) im*Ep*kx im*Ep*ky im*Ep*kz;
        -1*im*Ep*kx -1*L*kx^2-M*(ky^2+kz^2) -1*N*kx*ky -1*N*kx*kz;
        -1*im*Ep*ky -1*N*kx*ky -1*L*ky^2-M*(kx^2+kz^2) -1*N*ky*kz;
        -1*im*Ep*kz -1*N*kx*kz -1*N*ky*kz -1*L*kz^2-M*(kx^2+ky^2)]
    return h4
end
function soc(delta)
    hsoc=zeros(8,8)im;
    hsoc[2,3]=-1*im*(delta/3); hsoc[2,8]=delta/3; hsoc[3,2]=1*im*(delta/3); hsoc[3,8]=-1*im*(delta/3);
    hsoc[4,6]=-1*(delta/3); hsoc[4,7]=1*im*(delta/3); hsoc[6,4]=-1*(delta/3); hsoc[6,7]=1*im*(delta/3);
    hsoc[7,4]=-1*im*(delta/3); hsoc[7,6]=-1*im*(delta/3); hsoc[8,2]=(delta/3); hsoc[8,3]=1*im*(delta/3);
    return hsoc
end

function HsimpleSoc(Eg,kx,ky,kz,A,B,C,D,delta,F)
    htotS=zeros(8,8)im;
    h=Hsimple(Eg,kx,ky,kz,A,B,C,D);
    htotS[1:4,1:4]=h; htotS[5:8,5:8]=h;
    hsoc=soc(delta);
    htot=htotS+hsoc
    return htot
end
function HKaneSoc(kx,ky,kz,Eg,Ep,L,N,M,delta,F)
    htotS=zeros(8,8)im;
    h=HKane(kx,ky,kz,Eg,Ep,L,N,M,F);
    htotS[1:4,1:4]=h; htotS[5:8,5:8]=h;
    hsoc=soc(delta)
    htot=htotS+hsoc
    return htot
end 

