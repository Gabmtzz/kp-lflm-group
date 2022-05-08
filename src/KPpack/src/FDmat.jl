function mesh(Npts,layer,n)
    len=0
    boundary=zeros(n);
    
    for i in 1:length(layer)
        len += layer[i].size
        boundary[i]=len
    end
    dx=len/Npts; X=zeros(Npts)
    
    for i in 1:Npts
        X[i]=(i-1)*dx
    end
    
    return X, boundary

end

function createFuncFD(matriz,Emomentum)
    paramsFunc=["g_1","g_2", "g_3", "E_g","E_p", "F","k","Δ","VBO","c","cp"]; paramsFunSymb=functionParams(paramsFunc,Emomentum);
    funcMre=build_function(real(matriz),paramsFunSymb, expression=Val{false}); funcMim=build_function(imag(matriz),paramsFunSymb, expression=Val{false})
    mRe=funcMre[2]; mIm=funcMim[2];
    
    return mRe, mIm;
end

function evalFuncFD(mRe,mIm,mlayer,i,kx,ky,consth,const2)
    mq=mlayer[i]
    mReal=zeros(8,8); mImag=zeros(8,8)
    mRe(mReal,[kx,ky,0.0,mq.g1,mq.g2,mq.g3,mq.Eg,sqrt(mq.Ep),mq.F,mq.k,mq.delta,mq.VBO,consth,const2]); 
    mIm(mImag,[kx,ky,0.0,mq.g1,mq.g2,mq.g3,mq.Eg,sqrt(mq.Ep),mq.F,mq.k,mq.delta,mq.VBO,consth,const2])
    mEval=mReal+im*mImag;
    
    return mEval
end

function createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2)
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
            
    H0i=evalFuncFD(H0Re,H0Im,mlayer,i,kx,ky,consth,const2)
    H2im1=evalFuncFD(H2Re,H2Im,mlayer,im1,kx,ky,consth,const2); H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2); H2ip1=evalFuncFD(H2Re,H2Im,mlayer,ip1,kx,ky,consth,const2);
    
    Am=-(H2im1+2*H2i+H2ip1)/(2*(dx^2))+H0i; 
    
    return Am
end

function createFDBmatrix(mlayer,i,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
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
    
    H1i=evalFuncFD(H1Re,H1Im,mlayer,i,kx,ky,consth,const2); H1ip1=evalFuncFD(H1Re,H1Im,mlayer,ip1,kx,ky,consth,const2);
    H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2); H2ip1=evalFuncFD(H2Re,H2Im,mlayer,ip1,kx,ky,consth,const2);
    
    Bm=(H2ip1+H2i)/(2*(dx^2))+(H1ip1-H1i')/(2*dx);
    
    return Bm
end

function createFDCmatrix(mlayer,i,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
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
    
    H1im1=evalFuncFD(H1Re,H1Im,mlayer,im1,kx,ky,consth,const2); H1i=evalFuncFD(H1Re,H1Im,mlayer,i,kx,ky,consth,const2);
    H2im1=evalFuncFD(H2Re,H2Im,mlayer,im1,kx,ky,consth,const2); H2i=evalFuncFD(H2Re,H2Im,mlayer,i,kx,ky,consth,const2);
    
    Cm=(H2im1+H2i)/(2*(dx^2))-(H1im1-H1i')/(2*dx);
    
    return Cm
end

function FDHamiltonian(H0,H1,H2,mlayer,kx,ky,dx,consth,const2,len,Emomentum)
    H0Re,H0Im=createFuncFD(H0,Emomentum); H1Re,H1Im=createFuncFD(H1,Emomentum); H2Re,H2Im=createFuncFD(H2,Emomentum);
    
    AmV=Array{Matrix{ComplexF64}}(undef, len); CmV=Array{Matrix{ComplexF64}}(undef, len-1); BmV=Array{Matrix{ComplexF64}}(undef, len-1)  

    for i in 1:len
        if i==1
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2)
            BmV[i]=createFDBmatrix(mlayer,i+1,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
        elseif i==len
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2)
            CmV[i-1]=createFDCmatrix(mlayer,i-1,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
        else
            AmV[i]=createFDAmatrix(mlayer,i,H0Re,H0Im,H2Re,H2Im,kx,ky,dx,consth,const2)
            BmV[i]=createFDBmatrix(mlayer,i+1,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
            CmV[i-1]=createFDCmatrix(mlayer,i-1,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
        end
    end
    hns=Matrix(BlockTridiagonal(CmV, AmV, BmV))
    lhm=size(hns)[1]
    hns[1:8,lhm-7:lhm]=createFDCmatrix(mlayer,1,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
    hns[lhm-7:lhm,1:8]=createFDBmatrix(mlayer,len,H1Re,H1Im,H2Re,H2Im,kx,ky,dx,consth,const2)
    
    #hsp=sparse(hns)
    
    return sparse(hns)
end