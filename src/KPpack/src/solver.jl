function calcBandBulk(mm,kmax,Nt,model,Emomentum,constC,constCP)
    n=size(model)[1]

    #X
    pl=1; pm=0; pn=0
    Ex, Kx =KPpack.DiagM(mm,kmax+0.4,Nt,pl,pm,pn,model,Emomentum,n,constC,constCP);
    #L
    pl=1; pm=1; pn=1; Nt=100;
    El, Kl =KPpack.DiagM(mm,kmax,Nt,pl,pm,pn,model,Emomentum,n,constC,constCP);
    Exr=vcat(transpose(Kx),transpose(Ex))'
    el=vcat(transpose(-1*Kl),transpose(El))'
    Elr=rotl90(el,3)';
    dat=vcat(Elr,Exr);
    Ktot=dat[:,1]; Etot=dat[:,2:end];
    return Ktot,Etot
end

# ================================================================================
# function DiagM(mm,kmax,Nt,pl,pm,pn,soc,tipo)  
# m (structure Materials) => structure materials variable; kmax (float) => maximum value of k coordinate
# pl, pm, pn (integer) => numbers similar to Miller index [l,m,n]; soc (boolean) => enables spin-orbit calculation
# tipo (string) => specifies the Hamiltonian used in the calculation "Kane", "Simple"
#
# This function calculate the eigen values of the kp hamiltonian
# Return=> two arrays, one with the energies and other with the k values ===> dim 
# ================================================================================
function DiagM(mm,kmax,Nt,pl,pm,pn,model,Emomentum,n,constC,constCP)

    paramsFunc=["g_1","g_2", "g_3", "E_g","E_p", "F","k","Δ","VBO","c","cp","P"]
    
    paramsFunSymb=functionParams(paramsFunc,Emomentum)
    
    FuncMatreal=build_function(real(model),paramsFunSymb, expression=Val{false}); FuncMatim=build_function(imag(model),paramsFunSymb, expression=Val{false});
    HBulkRe=FuncMatreal[2]; HBulkIm=FuncMatim[2];

    En=zeros(Nt,n); Kp=zeros(Nt);
    for Nk in 1:Nt
       k=[pl,pm,pn]*kmax*(Nk-1)/(Nt+1); # create the k vector 
       kx=k[1]; ky=k[2];  kz=k[3]; #components of vector k

       HbRe=zeros(n,n); HbIm=zeros(n,n) 
       HBulkRe(HbRe,[kx,ky,kz,mm.g1,mm.g2,mm.g3,mm.Eg,mm.Ep,mm.F,mm.k,mm.delta,mm.VBO,constC,constCP,sqrt(mm.Ep)])
       HBulkIm(HbIm,[kx,ky,kz,mm.g1,mm.g2,mm.g3,mm.Eg,mm.Ep,mm.F,mm.k,mm.delta,mm.VBO,constC,constCP,sqrt(mm.Ep)])

       h=HbRe+im*HbIm


       w=eigvals(h); E=real(sort(w));
       En[Nk,:]=E;
       Kp[Nk]=(Nk-1)/(Nt+1);
    end
    return En,Kp
end

# =================================================================================
# function DOS(Ein,Eend,Estep,E,g)
# Ein (float)=>  inital energy;    Eend (float)=> final energy;   Estep (float)=> step energy for calculation
#
# Calculates a Density Of States between Ein and Eend with an energy separation of Estep
#
# Return => two arrays with the energy values and the Density of States values
# ================================================================================
function DOS(Ein,Eend,Estep,E,g)
    Edos=collect(range(Ein, step=Estep, stop=Eend));
    de=Edos[2]-Edos[1];
    ArrDos=zeros(length(Edos));
    for En=1:length(Edos)
        for i=1:size(E)[1]
            for j=1:size(E)[2]
                ArrDos[En]+=(g/(2*pi))/((Edos[En]-E[i,j])^2+(g/2)^2);
            end
        end
    end
    ArrDos=ArrDos/(de*sum(ArrDos))
    return Edos, ArrDos
end 

function DiagQWM(mlayer,kmax,Nt,dx,pl,pm,Npts,H0,H1r,H1l,H2,nc,nv,c,cp,sV,sC,Emomentum,pb,cr,s)
    En=zeros(Nt,nc); Env=zeros(Nt,nv); Kp=zeros(Nt);
    for Nk in 1:Nt
        k=[pl,pm]*kmax*(Nk-1)/(Nt+1);
        kx=k[1]; ky=k[2]
        Hamqw=FDHamiltonian(H0,H1r,H1l,H2,mlayer,kx,ky,dx,c,cp,Npts,Emomentum,pb,cr,s);
        λ, ϕ = eigs(Hamqw, nev=nc, which=:LM, sigma=sC);
        λ1, ϕ = eigs(Hamqw, nev=nv, which=:LM, sigma=sV);
        E=sort(real(λ));
        Ev=sort(real(λ1));
        En[Nk,:]=E;
        Env[Nk,:]=Ev;
        Kp[Nk]=(Nk-1)/(Nt+1)*kmax;
    end

    return En,Env,Kp
end

function EigSolQW(mlayer,Npts,H0,H1r,H1l,H2,c,cps,dx,Emomentum,pb,cr,s)
    kx,ky=0.0,0.0;
    Hamqw=FDHamiltonian(H0,H1r,H1l,H2,mlayer,kx,ky,dx,c,cps,Npts,Emomentum,pb,cr,s);
    Eqw0, EVqw0 = eigen(Matrix(Hamqw));
    return Eqw0, EVqw0
end

function eigenValQW(elem)
    inf,sup=0,0;
    pos=0
    for i in 1:length(elem)-1
        if (elem[i]<=0.0 && elem[i+1]>0.0)
            inf=elem[i]; sup=elem[i+1]
            pos=i
        end
    end
    return inf, sup, pos
end