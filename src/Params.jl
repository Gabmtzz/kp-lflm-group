# ================================================================
# Function params (Materials)
# Materials (Materials struct) => structure materials variable
# Return the kp parameters of a material from csv file.
# ================================================================

function params(Materials)
    arr=readdlm("MaterialsBulk.csv",',');
    index=findall(x->x==Materials.material,arr);
    i=index[1][1];
    Materials.g1=arr[i,2]; Materials.g2=arr[i,3]; Materials.g3=arr[i,4]; 
    Materials.Eg=arr[i,5]; Materials.Ep=arr[i,7]; Materials.F=arr[i,8]; 
    Materials.delta=arr[i,6];
    Materials.VBO=arr[i,9]; Materials.k=arr[i,10]; Materials.me=arr[i,11]
    Materials.al=arr[i,12];  Materials.B=arr[i,13]

    nothing 
end

# ============================================================================
# Function TempPar(material)
# Materials (Materials struct) => structure materials variable
# 
# obtains the Varshni parameter of respective material indicated in Materials 
# variable.
#
# return; two varables (float) with the values of the alpha and beta
# ============================================================================
function TempPar(material)
    arr1=readdlm("TempPar.csv",',');
    index1=findall(x->x==material,arr1);
    i=index1[1][1];
    alfa=arr1[i,2]; beta=arr1[i,3]; alat=arr1[i,4]
    
    return alfa, beta, alat
end

# ====================================================================================
# function aEgTemp(T,material,Eg)
# T (float)=> Temperature,  material (Materials struct) =>structure materials variable
# Eg (float) => energy gap
#
# Calculates the band gap energy for the temperature (T) using the Varshni formula
#
# Return EgT => Energy band gap 
# ====================================================================================

function aEgTemp(T,material,Eg,al)
    mat1,mat2,comp,alloy=DetMat(material)
    if mat1=="HgTe"
        EgT=-0.303
    else
        if mat2==""
            alfa, beta,alb = TempPar(mat1)
        else
            alfa1, beta1, alb1 = TempPar(mat1)
            alfa2, beta2, alb2 = TempPar(mat2)
        
            alfa=comp*alfa1+(1.0-comp)*alfa2
            beta=comp*beta1+(1.0-comp)*beta2
            alb=comp*alb1+(1.0-comp)*alb2
        end
        alT=al+alb*(T-300)
        EgT=Eg-(alfa*T^2)/(T+beta);
    end
    return EgT,alT
end

function getF(mm,opt)
    if opt=="calc"
        S=(1/mm.me)-(mm.Ep*(mm.Eg+(2/3)*mm.delta)/(mm.Eg*(mm.Eg+mm.delta)))

        mm.F=(S-1)/2
    elseif opt=="zero"
        S=0
        mm.F=(S-1)/2

        mm.Ep=(1/mm.me)*((mm.Eg*(mm.Eg+mm.delta))/(mm.Eg+(2/3)*mm.delta))
    elseif opt=="one"
        S=1.0
        mm.F=(S-1)/2

        mm.Ep=((1/mm.me)-1)*((mm.Eg*(mm.Eg+mm.delta))/(mm.Eg+(2/3)*mm.delta))
    else opt=="table"

    end



    return mm
end

# =======================================================================================
# function detMat(mat)
# mat (String) => Name of Alloy or binary compound
#
# Returns the the alloy's compositions and the respective concentration, if the alloy is 
# a binary compound this function returns only the name and the conmosition is zero.
#
# return: mat1, mat2 (string)=> name of binary compounds forming the alloy, 
# comp (float) => composition fraction, alloy (string) => name of the alloy
# =======================================================================================

function DetMat(mat)
    if length(mat)==4
        mat1=mat; mat2=""; comp=0.0; alloy="";
    else
        mat1="$(mat[1:2])$(mat[5:6])";mat2="$(mat[3:4])$(mat[5:6])"; comp=parse(Float64,mat[8:end]);alloy=mat[1:6];
    end
    return mat1, mat2, comp, alloy
end

# =======================================================================================
# function BowingPar (bowpar)
# bowpar (struct BowPar) => conatins the name of alloy and the bowing parameters
#
# This function adds the values of bowing parameters to BowPar structure
# =======================================================================================

function BowingPar(bowpar)
    arr=readdlm("BowingPar.csv",',');
    index1=findall(x->x==bowpar.alloy,arr);
    i=index1[1][1];
    bowpar.cEg1=arr[i,2]; bowpar.cEg2=arr[i,3]; bowpar.cEp=arr[i,5]; bowpar.cF=arr[i,6]; bowpar.cDelta=arr[i,4];
    bowpar.cVBO=arr[i,7]; bowpar.cme=arr[i,8] 
    nothing
end

# =======================================================================================
# function ParrAll(material1,material2,allMat,comp,All1)
# material1, material2 (string)=> name of binary compounds forming the alloy, 
# comp (float) => composition fraction, allMat (struct Materials) => struct with the data
# about the alloy
# All1(struct BowPar) => conatins the name of alloy and the bowing parameters
#
# this function calculates the parameters used in calculations for an alloy
# =======================================================================================

function ParrAll(material1,material2,allMat,comp,All1)
    allMat.g1=comp*material1.g1+(1.0-comp)*material2.g1; allMat.g2=comp*material1.g2+(1.0-comp)*material2.g2;
    allMat.g3=comp*material1.g3+(1.0-comp)*material2.g3;
    allMat.Eg=comp*material1.Eg+(1.0-comp)*material2.Eg-comp*(1-comp)*(All1.cEg1+All1.cEg2*comp);
    allMat.Ep=comp*material1.Ep+(1.0-comp)*material2.Ep-comp*(1-comp)*All1.cEp;
    allMat.F=comp*material1.F+(1.0-comp)*material2.F-comp*(1-comp)*All1.cF;
    allMat.delta=comp*material1.delta+(1.0-comp)*material2.delta-comp*(1-comp)*All1.cDelta;
    allMat.VBO=comp*material1.VBO+(1.0-comp)*material2.VBO-comp*(1-comp)*All1.cVBO;
    allMat.k=comp*material1.k+(1.0-comp)*material2.k;
    allMat.me=comp*material1.me+(1.0-comp)*material2.me-comp*(1-comp)*All1.cme
    allMat.al=comp*material1.al+(1.0-comp)*material2.al
    allMat.B=comp*material1.B+(1.0-comp)*material2.B
    nothing
end

# =======================================================================================
# function ParMat(AllMat,T)
# allMat (struct Materials) => struct with the data of the alloy
# T (float) => Temparature
# 
# This function fills the struct materials used in K.P bulk calculation
# =======================================================================================

function ParMat(AllMat,T,opt)
    mat1, mat2, comp, alloy = DetMat(AllMat.material)
    if mat2==""
        params(AllMat);
    else
        material1=Materials(mat1,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0); material2=Materials(mat2,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        params(material1); params(material2); 
        all1=BowPar(alloy,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        BowingPar(all1);
        ParrAll(material1,material2,AllMat,comp,all1)
    end
    EgT,alT=aEgTemp(T,AllMat.material,AllMat.Eg,AllMat.al)
    AllMat.Eg=EgT; AllMat.al=alT
    AllMat=getF(AllMat,opt)
end

# =======================================================================================
# function supParams(layer,X,boundary,mlayer,T)
# 
# layer ()
# =======================================================================================

function supParams(layer,X,boundary,mlayer,T,opt)
    nlay=1;
    boundaryPoints=zeros(length(boundary));
    for i in 1:length(X)
        mlayer[i]=Materials(layer[nlay].material,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        ParMat(mlayer[i],T,opt)
        mlayer[i].Eg=mlayer[i].Eg
        if X[i]>= boundary[nlay] boundaryPoints[nlay]=i; nlay+=1  end
    end    
    return boundaryPoints[1:length(boundary)-1]
end

function setStructure(structure,fB,Temp)
    nlayer=size(structure)[1]
    alLay=Array{Float64}(undef,nlayer)
    dis=parse.(Float64,structure[:,2]); xdisTot=sum(dis)

    layer=Array{mat}(undef,nlayer,)

    for i in 1:nlayer
        layer[i]=mat(structure[i,1],dis[i]);
        mml=Materials(structure[i,1],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        ParMat(mml,Temp,"table");
        alLay[i]=(dis[i]/xdisTot)*mml.al
    end
    aprom=sum(alLay)
    Npts=Int(floor(2*fB*(xdisTot/aprom)))

    return Npts+1,layer,nlayer
end


function writefile(data,name)
    path="../../data/$(name).csv"
    fout=open(path,"w")
    writedlm(fout, data, ",", header=false)
    close(fout)
end

# =================================================================================================================
# The data come from:
# * Igor Vurgaftman, Matthew P. Lumb, and Jerry R. Meyer, Bands and Photons in III-V Semiconductor Quantum Structures,Oxford, 2020
# 
# *Michał Marchewka, et.al., Finite-difference method applied for eight-band kp model for Hg1−xCdxTe/HgTe quantum well, International Journal of 
#  Modern Physics B, 2017
#
# * I. Vurgaftman and J. R. Meyer , Band parameters for III–V compound semiconductors and their alloys", Journal of Applied Physics 89, 5815-5875 # (2001)
# =================================================================================================================

