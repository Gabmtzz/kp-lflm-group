function PlotBands(Etot,Ktot)
    plt.figure()
    plt.xticks(ticks=[-1.0,-0.5,0,0.5,1.0], labels=[L"$\leftarrow$ L","Λ"," Γ", "Δ",L"X $\rightarrow$"])
    plt.ylabel("Energy [eV]")
    plt.plot(Ktot,Etot, color="black")
end

function PlotDOS(EDOS,aDOS)
    plt.figure()
    plt.yticks([])
    plt.xlabel("Energy [eV]"); plt.ylabel("DOS")
    plt.plot(EDOS,aDOS, color="black")
end

function plotProf(mlayer,X,option)
    Ec=zeros(length(mlayer)); Ev=zeros(length(mlayer)) 
    for i in 1:length(mlayer)
        Ec[i]=mlayer[i].Eg+mlayer[i].VBO
        Ev[i]=mlayer[i].VBO
    end
    
    plt.xlabel("X [nm]"); plt.ylabel("Energy [Ev]")
    if option=="both"
        plt.plot(X,Ec, color="blue", label="Ec")
        plt.plot(X,Ev, color="red", label="Ev")
        plt.legend()
    elseif option=="Ec"
        plt.plot(X,Ec, color="black")
    elseif option=="Ev"
        plt.plot(X,Ev, color="black")
    else 
        print("bad option: use both, Ec or Ev")
    end
end

function PlotQWBand(Ecq11,Evq11,Kqw11,Ecq10,Evq10,Kqw10,option,poslab, kmax)
    plt.xlabel(L"$k_{||}~ [nm^{-1}]$"); plt.ylabel("Energy [Ev]")
    if option=="complete_Band"
        plt.text(-0.7*kmax,poslab, L"$\leftarrow~~[110]$"); plt.text(0.7*kmax,poslab, L"$[100]~~\rightarrow$")
        plt.plot(-1*Kqw11,Ecq11,color="black" )
        plt.plot(-1*Kqw11,Evq11,color="black" )
        plt.plot(Kqw10,Evq10,color="black" )
        plt.plot(Kqw10,Ecq10,color="black" )
    elseif option=="complete_C"
        plt.text(-0.7*kmax,poslab, L"$\leftarrow~~[110]$"); plt.text(0.7*kmax,poslab, L"$[100]~~\rightarrow$")
        plt.plot(-1*Kqw11,Ecq11,color="black" )
        plt.plot(Kqw10,Ecq10,color="black" )
    elseif option=="complete_V"
        plt.text(-0.7*kmax,poslab, L"$\leftarrow~~[110]$"); plt.text(0.7*kmax,poslab, L"$[100]~~\rightarrow$")
        plt.plot(-1*Kqw11,Evq11,color="black" )
        plt.plot(Kqw10,Evq10,color="black" )
    elseif option=="comp_V"
        plt.plot(Kqw11,Evq11,color="black", label=L"$[110]$")
        plt.plot(Kqw10,Evq10,color="red",linestyle="dotted", label=L"$[100]$")
        #plt.legend()
    elseif option=="comp_C"
        plt.plot(Kqw11,Ecq11,color="black" , label=L"$[110]$")
        plt.plot(Kqw10,Ecq10,color="red",linestyle="dotted", label=L"$[100]$" )
        #plt.legend()
    else 
        print("bad option: use complete_Band, complete_C, complete_V, comp_V or comp_C")
    end
end

function PloteigvQW(Npts,Eqw0,siz)
    x1=collect(range(1,siz*Npts, length=siz*Npts));
    Egqw=real(Eqw0);
    plt.xlabel("Eigenvalue number, α"); plt.ylabel("Energy [eV]")
    plt.plot(x1,Egqw, color="blue", "o")
end

function plotEigSurface(EArr,eig,kinic,kend)
    plt.ylabel(L"k_x ~[nm^{-1}]"); plt.xlabel(L"k_y~[nm^{-1}]") 
    plt.xticks(ticks=collect(range(0,step=2.5,stop=20)),labels=collect(range(kinic,step=2/8,stop=kend)))
    plt.yticks(ticks=collect(range(0,step=2.5,stop=20)),labels=collect(range(kinic,step=2/8,stop=kend)))
    plot=plt.imshow(EArr[:,:,eig],cmap="turbo",aspect="auto",interpolation="gaussian")
    bar=plt.colorbar(plot)
    bar.set_label(L"Energy~[eV]")
end

function plotProbDistr(iVec,pos,EVqw0,siz,vecband,boundary,X,legendArr)
    plt.xlabel("X [nm]"); plt.ylabel(L"|\Psi|^2 ~~[U.A.]")
    vv,Evec=zeros(length(X)),zeros(length(X))*im
    valMaxArr1,valMaxArr2=zeros(length(iVec)),zeros(length(vecband))
    for k in 1:length(vecband)
        for i in 1:length(iVec)
            Evec=EVqw0[vecband[k]:siz:end,pos+iVec[i]]

            vv=abs.(Evec)
            plt.plot(X,real(vv),label=legendArr[k])
            valMaxArr1[i]=sort(vv,rev=true)[1]
        end
        valMaxArr2[k]=sort(valMaxArr1,rev=true)[1]
        
    end
    plt.legend()
    dx=X[2]-X[1]
    valmax=sort(valMaxArr2,rev=true)[1]
    for j in 1:length(boundary)-1
        l0i,l0f=[boundary[j]+dx,boundary[j]+dx],[0.0,valmax]; plt.plot(l0i,l0f,color="black", linestyle="--")
    end
end