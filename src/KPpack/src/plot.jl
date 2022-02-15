ENV["GKS_ENCODING"]="utf-8" 
using Plots

function PlotBands(Ex,El,Kx,Kl)
    plot(Kx,Ex,color=:blue, leg=false)
    plot!(-1*Kl,El,color=:blue, leg=false)
    plot!(ylabel="Energy (EV)", xticks=([-1.0,-0.5,0,0.5,1.0],["<- L","Λ","Γ","Δ","X ->"]))
end
function PlotDOS(EDOS,aDOS)
    plot(EDOS,aDOS,color=:blue, leg=false)
    plot!(xlabel="Energy (Ev)", yticks=([],[]),ylabel="DOS")
end

function plotProf(mlayer,X,option)
    Ec=zeros(length(mlayer)); Ev=zeros(length(mlayer)) 
    for i in 1:length(mlayer)
        Ec[i]=mlayer[i].Eg
        Ev[i]=mlayer[i].VBO
    end
    if option=="both"
        plot(X,Ec, color=:blue, label="Ec")
        plot!(X,Ev, color=:red, label="Ev")
    elseif option=="Ec"
        plot(X,Ec, color=:blue, leg=false)
    elseif option=="Ev"
        plot(X,Ev, color=:blue, leg=false)
    end
    plot!(xlabel="X (nm)", ylabel="Energy (Ev)")
end

function PlotQWBand(Ecq11,Evq11,Kqw11,Ecq10,Evq10,Kqw10,option,yinf,ysup)
    if option=="complete_Band"
        plot(-1*Kqw11,Ecq11,color=:blue, leg=false)
        plot!(-1*Kqw11,Evq11,color=:blue, leg=false)
        plot!(Kqw10,Evq10,color=:blue, leg=false)
        plot!(Kqw10,Ecq10,color=:blue, leg=false)
    elseif option=="complete_C"
        plot(-1*Kqw11,Ecq11,color=:blue, leg=false)
        plot!(Kqw10,Ecq10,color=:blue, leg=false)
    elseif option=="complete_V"
        plot(-1*Kqw11,Evq11,color=:blue, leg=false)
        plot!(Kqw10,Evq10,color=:blue, leg=false)
    elseif option=="comp_V"
        plot(Kqw11,Evq11,color=:blue, leg=false)
        plot!(Kqw10,Evq10,color=:red,line=(:dot, 4), leg=false)
    elseif option=="comp_C"
        plot(Kqw11,Ecq11,color=:blue, leg=false)
        plot!(Kqw10,Ecq10,color=:red,line=(:dot, 4), leg=false)
    end
    plot!(ylabel="Energy (Ev)", xlabel="k (1/nm)", ylims=(yinf,ysup))
end