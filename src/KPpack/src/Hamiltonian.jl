function Hsimple(Eg,kx,ky,kz,A,B,C,D)
    h=[Eg+D*(kx^2+ky^2+kz^2) 0 0 0 ;
        0 -A*kx^2-B*(ky^2+kz^2) -C*kx*ky -C*kx*kz ;
        0 -C*kx*ky -A*ky^2-B*(kx^2+kz^2) -C*ky*kz ;
        0 -C*kx*kz -C*ky*kz -A*kz^2-B*(kx^2+ky^2)]
    return h
end
function HKane(kx,ky,kz,Eg,Ep,L,N,M,F,VBO)
    h4=[Eg+VBO+(1+2*F)*(kx^2+ky^2+kz^2) im*Ep*kx im*Ep*ky im*Ep*kz;
        -1*im*Ep*kx VBO-1*L*kx^2-M*(ky^2+kz^2) -1*N*kx*ky -1*N*kx*kz;
        -1*im*Ep*ky -1*N*kx*ky VBO-1*L*ky^2-M*(kx^2+kz^2) -1*N*ky*kz;
        -1*im*Ep*kz -1*N*kx*kz -1*N*ky*kz VBO-1*L*kz^2-M*(kx^2+ky^2)]
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
function HKaneSoc(kx,ky,kz,Eg,Ep,L,N,M,delta,F,VBO)
    htotS=zeros(8,8)im;
    h=HKane(kx,ky,kz,Eg,Ep,L,N,M,F,VBO);
    htotS[1:4,1:4]=h; htotS[5:8,5:8]=h;
    hsoc=soc(delta)
    htot=htotS+hsoc
    return htot
end 
function HKaneII(kx,ky,kz,m)
   kpar2=kx*kx+ky*ky; kplus=kx+im*ky; kminus=kx-im*ky; c=1; P=sqrt(c*m.Ep);
   mu=(m.g3-m.g2)/2; gamma=(m.g3+m.g2)/2 
   T=m.Eg+m.VBO+c*((2*m.F+1)*kpar2+(2*m.F+1)*kz*kz); U=m.VBO-c*m.g1*(kpar2+kz*kz)
   V=-1*c*m.g2*(kpar2-2*kz*kz); R=-c*sqrt(3)*(mu*kplus*kplus-gamma*kminus*kminus)
   Splus=-2*c*sqrt(3)*m.g3*kplus*kz; Sminus=-2*c*sqrt(3)*m.g3*kminus*kz

   hvec=[T,T,U+V,U-V,U-V,U+V,U-m.delta,U-m.delta]; hdiag=Diagonal(hvec)
    
   h=zeros(8,8)*im
   h[1,3]=-(1/sqrt(2))*P*kplus; h[1,4]=(sqrt(2)/sqrt(3))*P*kz; h[1,5]=(1/sqrt(6))*P*kminus 
   h[1,7]=-(1/sqrt(3))*P*kz; h[1,8]=-(1/sqrt(3))*P*kminus
    
   h[2,4]=-(1/sqrt(6))*P*kplus; h[2,5]=(sqrt(2)/sqrt(3))*P*kz; h[2,6]=(1/sqrt(2))*P*kminus
   h[2,7]=-(1/sqrt(3))*P*kplus; h[2,8]=(1/sqrt(3))*P*kz
    
   h[3,4]=-Sminus; h[3,5]=R; h[3,7]=(1/sqrt(2))*Sminus; h[3,8]=-sqrt(2)*R
    
   h[4,6]=R; h[4,7]=sqrt(2)*V; h[4,8]=-sqrt(3/2)*Sminus
    
   h[5,6]=Splus'; h[5,7]=-sqrt(3/2)*Splus; h[5,8]=-sqrt(2)*V
    
   h[6,7]=sqrt(2)*R'; h[6,8]=(1/sqrt(2))*Splus
    
   Htot=hdiag+h+h'
    
   return Htot 
end

# ===================================================================================================
# The hamiltonians are taken from:
# * Igor Vurgaftman, Matthew P. Lumb, and Jerry R. Meyer, Bands and Photons in III-V Semiconductor Quantum Structures,Oxford, 2020
# 
# *Michał Marchewka, et.al., Finite-difference method applied for eight-band kp model for Hg1−xCdxTe/HgTe quantum well, International Journal of 
#  Modern Physics B, 2017
# ===================================================================================================