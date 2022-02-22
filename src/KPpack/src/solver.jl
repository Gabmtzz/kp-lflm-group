# ==========================================================================================================================
# function DiagM(mm,kmax,Nt,pl,pm,pn,soc,tipo)  
# m (structure Materials) => structure materials variable; kmax (float) => maximum value of k coordinate
# pl, pm, pn (integer) => numbers similar to Miller index [l,m,n]; soc (boolean) => enables spin-orbit calculation
# tipo (string) => specifies the Hamiltonian used in the calculation "Kane", "Simple"
#
# This function calculate the eigen values of the kp hamiltonian
# Return=> two arrays, one with the energies and other with the k values ===> dim 
# ==========================================================================================================================
function DiagM(mm,kmax,Nt,pl,pm,pn,soc,tipo)
    n=1;
    L=mm.g1+4*mm.g2-(mm.Ep/mm.Eg); M=mm.g1-2*mm.g2; N=6*mm.g3-(mm.Ep/mm.Eg);
    A=L+(mm.Ep/mm.Eg); B=M; C=N+(mm.Ep/mm.Eg); D=(1+2*mm.F+(mm.Ep/mm.Eg));   #calculate the k.p parameters
    if soc n=2 end
    En=zeros(Nt,4*n); Kp=zeros(Nt);
    for Nk in 1:Nt
       k=[pl,pm,pn]*kmax*(Nk-1)/(Nt+1); # create the k vector 
       kx=k[1]; ky=k[2];  kz=k[3]; #components of vector k
       if soc  # verifies if the spin-orbit coupling is enabled
            if tipo=="Kane"
                
               h=HKaneSoc(kx,ky,kz,mm.Eg,sqrt(mm.Ep),L,N,M,mm.delta,mm.F,mm.VBO) # call the function HKaneSoc
            elseif tipo=="Simple"
               h=HsimpleSoc(mm.Eg,kx,ky,kz,A,B,C,D,mm.delta,mm.F)
            elseif tipo=="KaneII"
                 h=HKaneII(kx,ky,kz,mm)
                
            else
                
            end
                
        else
           if tipo=="Kane"
                
               h=HKane(kx,ky,kz,mm.Eg,sqrt(mm.Ep),L,N,M,mm.F,mm.VBO)
            elseif tipo=="Simple"
                h=Hsimple(mm.Eg,kx,ky,kz,A,B,C,D)
            elseif tipo=="KaneII"
                print("Not implemented")
                
            else
                
            end
        end
        w=eigvals(h); E=real(sort(w));
        En[Nk,:]=E;
        Kp[Nk]=(Nk-1)/(Nt+1);
    end
    return En,Kp
end

# ====================================================================================================
# function DOS(Ein,Eend,Estep,E,g)
# Ein (float)=>  inital energy;    Eend (float)=> final energy;   Estep (float)=> step energy for calculation
#
# Calculates a Density Of States between Ein and Eend with an energy separation of Estep
#
# Return => two arrays with the energy values and the Density of States values
# ====================================================================================================
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

function DiagQWM(mlayer,kmax,Nt,dx,pl,pm,Npts,boundPoints,n,c,cp,sV,sC )
    En=zeros(Nt,n); Env=zeros(Nt,n); Kp=zeros(Nt);
    for Nk in 1:Nt
        k=[pl,pm]*kmax*(Nk-1)/(Nt+1);
        kx=k[1]; ky=k[2]
        Hamqw=QWHamiltonianMatrix(mlayer,kx,ky,dx,Npts,boundPoints,c,cp);
        λ, ϕ = eigs(Hamqw, nev=n, which=:LM, sigma=sC);
        λ1, ϕ = eigs(Hamqw, nev=n, which=:LM, sigma=sV);
        E=sort(real(λ));
        Ev=sort(real(λ1));
        En[Nk,:]=E;
        Env[Nk,:]=Ev;
        Kp[Nk]=(Nk-1)/(Nt+1)*kmax;
    end

    return En,Env,Kp
end

function EigSolQW(mlayer,Npts,boundPoints,c,cps,dx)
    kx,ky=0.0,0.0;
    Hamqw=KPpack.QWHamiltonianMatrix(mlayer,kx,ky,dx,Npts,boundPoints,c,cps);
    Eqw0, EVqw0 = eigen(Matrix(Hamqw));
    return Eqw0, EVqw0
end