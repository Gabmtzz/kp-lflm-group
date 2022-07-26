function mesh(Npts,layer,n)
    len=0
    boundary=zeros(n);
    
    for i in 1:length(layer)
        len += layer[i].size
        boundary[i]=len
    end
    dx=len/Npts; X=zeros(Npts)
    
    for i in 1:Npts
        X[i]=(i)*dx
    end
    
    return X, boundary

end

function createFuncFD(matriz,Emomentum)
    paramsFunc=["g_1","g_2", "g_3", "E_g","E_p", "F","k","Δ","VBO","c","cp","P","s"]; paramsFunSymb=functionParams(paramsFunc,Emomentum);
    funcMre=build_function(real(matriz),paramsFunSymb, expression=Val{false}); funcMim=build_function(imag(matriz),paramsFunSymb, expression=Val{false})
    mRe=funcMre[2]; mIm=funcMim[2];
    
    return mRe, mIm;
end

function evalFuncFD(mRe,mIm,mlayer,i,kx,ky,consth,const2,cr,s,siz)
    mq=mlayer[i]
    mReal=zeros(siz,siz); mImag=zeros(siz,siz)
    mRe(mReal,[kx,ky,0.0,mq.g1,mq.g2,mq.g3,mq.Eg,cr*mq.Ep,mq.F,mq.k,mq.delta,mq.VBO,consth,const2,sqrt(cr*mq.Ep),s]); 
    mIm(mImag,[kx,ky,0.0,mq.g1,mq.g2,mq.g3,mq.Eg,cr*mq.Ep,mq.F,mq.k,mq.delta,mq.VBO,consth,const2,sqrt(cr*mq.Ep),s])
    mEval=mReal+im*mImag;
    
    return mEval
end

function createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
    if i==1
        im1=i
        ip1=i+1
    elseif i==length(mlayer)
        im1=i-1
        ip1=i
    else
        im1=i-1
        ip1=i+1
    end
    H0i=evalFuncFD(H0Re,H0Im,mlayer,i,kx,ky,consth,const2,cr,s,siz)
    H2im1=evalFuncFD(H2Re,H2Im,mlayer,im1,kx,ky,consth,const2,cr,s,siz); H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2,cr,s,siz); H2ip1=evalFuncFD(H2Re,H2Im,mlayer,ip1,kx,ky,consth,const2,cr,s,siz);
    
    Am=-(H2im1+2*H2i+H2ip1)/(2*(dx^2))+H0i; 
    
    return Am
end

function createFDBmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
    if i==1
        im1=i
        ip1=i+1
    elseif i==length(mlayer)
        im1=i-1
        ip1=i
    else
        im1=i-1
        ip1=i+1
    end
    
    H1i=evalFuncFD(H1lRe,H1lIm,mlayer,i,kx,ky,consth,const2,cr,s,siz); H1ip1=evalFuncFD(H1rRe,H1rIm,mlayer,ip1,kx,ky,consth,const2,cr,s,siz);
    H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2,cr,s,siz); H2ip1=evalFuncFD(H2Re,H2Im,mlayer,ip1,kx,ky,consth,const2,cr,s,siz);
    
    Bm=(H2ip1+H2i)/(2*(dx^2))+(H1ip1+H1i)/(2*dx);
    
    return Bm
end

function createFDCmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
    if i==1
        im1=i
        ip1=i+1
    elseif i==length(mlayer)
        im1=i-1
        ip1=i
    else
        im1=i-1
        ip1=i+1
    end
    H1im1=evalFuncFD(H1rRe,H1rIm,mlayer,im1,kx,ky,consth,const2,cr,s,siz); H1i=evalFuncFD(H1lRe,H1lIm,mlayer,i,kx,ky,consth,const2,cr,s,siz);
    H2im1=evalFuncFD(H2Re,H2Im,mlayer,im1,kx,ky,consth,const2,cr,s,siz); H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2,cr,s,siz);
    
    Cm=(H2im1+H2i)/(2*(dx^2))-(H1im1+H1i)/(2*dx);
    
    return Cm
end

function FDHamiltonian(H0,H1r,H1l,H2,mlayer,kx,ky,dx,consth,const2,len,Emomentum,pb,cr,s)
    siz=size(H0)[1] 
    H0Re,H0Im=createFuncFD(H0,Emomentum); H1rRe,H1rIm=createFuncFD(H1r,Emomentum); H1lRe,H1lIm=createFuncFD(H1l,Emomentum)
    H2Re,H2Im=createFuncFD(H2,Emomentum);
    
    AmV=Array{Matrix{ComplexF64}}(undef, len); CmV=Array{Matrix{ComplexF64}}(undef, len-1); BmV=Array{Matrix{ComplexF64}}(undef, len-1)  

    for i in 1:len
        if i==1
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
            BmV[i]=createFDBmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
        elseif i==len
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
            CmV[i-1]=createFDCmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
        else
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
            BmV[i]=createFDBmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
            CmV[i-1]=createFDCmatrix(mlayer,i,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
        end
    end
    hns=Matrix(BlockTridiagonal(CmV, AmV, BmV))
    if pb
        lhm=size(hns)[1]
        hns[1:siz,lhm-siz-1:lhm]=createFDCmatrix(mlayer,1,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
        hns[lhm-siz-1:lhm,1:siz]=createFDBmatrix(mlayer,len,H1rRe,H1rIm,H1lRe,H1lIm,H2Re,H2Im,kx,ky,dx,consth,const2,cr,s,siz)
    end
    
    #hsp=sparse(hns)
    
    return sparse(hns)
end