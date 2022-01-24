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