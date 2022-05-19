function substituteValue(symexp,varSub,value)
    realPart=real(symexp); imagPart=imag(symexp)
    subRe=substitute(realPart,Dict(varSub=>value)); subIm=substitute(imagPart,Dict(varSub=>value))
    subRe=simplify.(subRe); subIm=simplify.(subIm) 
    subVal=subRe+im*subIm
    
    return subVal
end

function extractParams(M)
    s=Symbolics.get_variables(M[1,1])
    for i in 1:size(M)[1]
        for j in 1:size(M)[2]
           if i==1 && j==1
                    s=Symbolics.get_variables(M[i,j])
                else
                    sAux=Symbolics.get_variables(M[i,j]);
                    s=union(s,setdiff(sAux,s))
                end 
        end
    end
    
    return s
end

function getVar(H)
    parHr=extractParams(real(H)); parHi=extractParams(imag(H));
    parH=union(parHr,parHi)

    se=Symbol("se"); se=@variables $se[1:length(parH)]; se=se[1]
    se=Symbolics.scalarize(se)

    for i in 1:length(se)
        se=substitute(se,Dict(se[i]=>parH[i]))
    end

    return se
end

function extractmatrix(matriz,var)
    Hz=substitute(matriz,Dict(var=>0))
    Hnk=matriz-Hz; Hnk=simplify.(Hnk)

    k0=Symbol("k0"); k0= @variables $k0[1:size(matriz)[1],1:size(matriz)[2]]; k0=k0[1] 
    k1=Symbol("k1"); k1= @variables $k1[1:size(matriz)[1],1:size(matriz)[2]]; k1=k1[1] 
    k2=Symbol("k2"); k2= @variables $k2[1:size(matriz)[1],1:size(matriz)[2]]; k2=k2[1] 

    k0=Symbolics.scalarize(k0); k1=Symbolics.scalarize(k1); k2=Symbolics.scalarize(k2) 


    for i in 1:size(matriz)[1]
        for j in 1:size(matriz)[2]
            v1=substituteValue(Hnk[i,j],var,-1); v2=substituteValue(Hnk[i,j],var,1)
            subst=simplify(v1+v2) 
            if typeof(subst==0)==typeof(true)
                if subst==0
                    k1[i,j]=substitute(k1[i,j],Dict(k1[i,j]=>Hnk[i,j]))
                    k2[i,j]=substitute(k2[i,j],Dict(k2[i,j]=>0))
                else
                    k2[i,j]=substitute(k2[i,j],Dict(k2[i,j]=>Hnk[i,j]))
                    k1[i,j]=substitute(k1[i,j],Dict(k1[i,j]=>0))
                end
            else
                k2[i,j]=substitute(k2[i,j],Dict(k2[i,j]=>Hnk[i,j]))
                k1[i,j]=substitute(k1[i,j],Dict(k1[i,j]=>0))
            end
        end
    end

    k0=simplify.(matriz-k1-k2)
    k1=substitute(k1,Dict(var => 1)); k2=substitute(k2,Dict(var => 1)); k0=substitute(k0,Dict(var => 1))
    
    return k0,k1,k2
end

function createKm(H,var)
    matrizr=real(H); matrizi=imag(H);
    k0r,k1r,k2r=extractmatrix(matrizr,var); k0i,k1i,k2i=extractmatrix(matrizi,var);
    k0=simplify.(k0r)+im*simplify.(k0i); k1=-im*simplify.(k1r)+simplify.(k1i); k2=-simplify.(k2r)-im*simplify.(k2i); 
    
    return k0,k1,k2
end

function constructMatrixFD(H0,H1,H2)
    Δz=Symbol("Δz"); Δz=@variables $Δz; Δz=Δz[1]

    A=-(4*H2/(2*Δz^2))+H0; B=(2*H2/(2*Δz^2))+((H1-H1')/(2*Δz)); C=(2*H2/(2*Δz^2))-((H1-H1')/(2*Δz));
    
    return A,B,C
end

function Paramsboundary(params,b)
    parBound=params
    for i in 1:length(parBound)
        symtS=params[i]
        Ess=Symbolics.tosymbol(symtS)
        nSyn=Symbol(Ess,"_",b)
        SymSubs= @variables $nSyn
        parBound[i]=SymSubs[1]
    end
    return parBound
end

function modifiesBCmatrix(matriz,PrH,exclude,b)
    params=setdiff(PrH,exclude)
    paramsB1=Paramsboundary(params,b)
    params=setdiff(PrH,exclude)

    matriz1=matriz

    for i in 1:length(paramsB1) 
        matriz1=substitute(matriz1,Dict(params[i]=>paramsB1[i]))
    end
    return matriz1
end
function createBoundarymatrix(matrix,PrH,exclude,b)
    matRe=real(matrix); matIm=imag(matrix)
    matReMod=modifiesBCmatrix(matRe,PrH,exclude,b); matImMod=modifiesBCmatrix(matIm,PrH,exclude,b)
    matMod = matReMod+im*matImMod
    
    return matMod 
end   

function createBmatrizFD(H0,H1,H2,PrH,exclude,b1,b2)
    
    Δz=Symbol("Δz"); Δz=@variables $Δz; Δz=Δz[1]
    H0l=createBoundarymatrix(H0,PrH,exclude,b2); H2l=createBoundarymatrix(H2,PrH,exclude,b2);  H2r=createBoundarymatrix(H2,PrH,exclude,b1);
    H0r=createBoundarymatrix(H0,PrH,exclude,b1); H1l=createBoundarymatrix(H1,PrH,exclude,b2);  H1r=createBoundarymatrix(H1,PrH,exclude,b1);
    
    Al=(-(H2r+3H2l)/(2Δz^2))+H0l; Ar=(-(H2l+3H2r)/(2Δz^2))+H0r
    Bl=((H2r+H2l)/(2Δz^2))+((H1r-H1l')/(2Δz)); Br=(2*H2r/(2*Δz^2))+((H1r-H1r')/(2*Δz))
    Cr=((H2l+H2r)/(2Δz^2))+((H1l-H1r')/(2Δz)); Cl=(2*H2l/(2*Δz^2))+((H1l-H1l')/(2*Δz))
    return Al,Ar,Bl,Br,Cr,Cl
end

function StrtoSymBase(expr)
    
    if typeof(expr)==Symbol
        ss=expr; arrv=@variables $ss; arrv=arrv[1]
        ex= Expr(:(=),expr,arrv); eval(ex)
    elseif typeof(expr)==Expr
        elems=expr.args
        inic=1
        if length(elems) ==2 
            inic=2
        end

        for i in inic:length(elems)
            el=elems[i]
            if el != :+ && el != :/&&el != :*&&el != :^ && el != :- && el != :sqrt 

                if typeof(el)==Symbol
                    ss=elems[i]; arrv=@variables $ss; arrv=arrv[1]
                    ex= Expr(:(=),elems[i],arrv); eval(ex)
                elseif typeof(el)==Expr
                    StrtoSymBase(el)
                end
            end
        end
    end
end

function StrtoSymbConv(s)
    sM=Meta.parse(s)


 #   for i in 1:length(sMArg)
 #       if sMArg[i] != :sqrt
            StrtoSymBase(sM)
  #      end
   # end
    strconv=eval(sM)
    
    return strconv
end

function StrtoSymbComplexConv(sc)
    strRe=""; strIm=""

    pos=findall.( "im", sc)

    if length(pos)!=0
        el=collect(pos[1])
    
        strRe=sc[1:el[1]-2]; strIm=sc[el[2]+1:end]
        exprRe=KPpack.StrtoSymbConv(strRe); exprIm=KPpack.StrtoSymbConv(strIm);
        expr=exprRe+im*exprIm
    else
        strRe=sc;
        expr=KPpack.StrtoSymbConv(strRe);
    end
        return expr
end

function createSymbMatrix(model)
    Hs=Symbol("Hs"); Hs= @variables $Hs[1:size(model)[1],1:size(model)[2]]; Hs=Hs[1]
    Hs=Symbolics.scalarize(Hs)

    HsIm=Symbol("HsIm"); HsIm= @variables $HsIm[1:size(model)[1],1:size(model)[2]]; HsIm=HsIm[1]
    HsIm=Symbolics.scalarize(HsIm)

    for i in 1:size(model)[1]
        for j in 1:size(model)[2]
            epxrS=model[i,j]; exSimb=KPpack.StrtoSymbComplexConv(epxrS)
            Hs[i,j]=substitute(Hs[i,j],Dict(Hs[i,j]=>real(exSimb))); HsIm[i,j]=substitute(HsIm[i,j],Dict(HsIm[i,j]=>imag(exSimb))) 
        end
    end
    Ht=Hs+im*HsIm
    return Ht;
end

function setMomentum(strK)
    k1=Symbol(strK[1]); k2=Symbol(strK[2]); k3=Symbol(strK[3]);

    momentumK = @variables $k1 $k2 $k3
    
    return momentumK
end

function functionParams(paramsFunc,momentum)
    parmsfunSym=Symbol("PsF"); parmsfunSym= @variables $parmsfunSym[1:length(paramsFunc)]; parmsfunSym=parmsfunSym[1] 

    parmsfunSym=Symbolics.scalarize(parmsfunSym)
    for i in 1:length(paramsFunc)
        simbF=Symbol(paramsFunc[i]); simbF = @variables $simbF; simbF=simbF[1]
        parmsfunSym[i]=substitute(parmsfunSym[i],Dict(parmsfunSym[i]=> simbF))
    end

    funcParams=union(momentum,parmsfunSym)
    return funcParams
end

function generateFunctionMatrix(matrix,paramsFunc,Emomentum)
    #path="./src/HamiltonianFunc/"
    
    paramsFunSymb=KPpack.functionParams(paramsFunc,Emomentum)

    FuncMatreal=build_function(real(matrix),paramsFunSymb); FuncMatim=build_function(imag(matrix),paramsFunSymb);
    
    #nameReal=path*name*"Real.jl"; nameImag=path*name*"Imag.jl" 
    
    #write(nameReal, string(FuncMatreal[2])); write(nameImag, string(FuncMatim[2]));
    return FuncMatreal[2], FuncMatim[2]
end

function getSfunct(name)
    path="./src/HamiltonianFunc/"
    nameReal=path*name*"Real.jl"; nameImag=path*name*"Imag.jl"
    
    fRe=open(nameReal,"r")
    sRe=read(fRe,String);
    close(fRe)
    
    fIm=open(nameImag,"r")
    sIm=read(fIm,String);
    close(fIm)
    
    return sRe, sIm
end

function setExcludePar(exclude)
    lEx=length(exclude)

    excludeS=Symbol("exS");excludeS= @variables $excludeS[1:lEx]; excludeS=excludeS[1]
    excludeS=Symbolics.scalarize(excludeS)

    for i in 1:lEx
        subs=Symbol(exclude[i]); subs= @variables $subs; subs=subs[1] 
        excludeS[i]=substitute(excludeS[i],Dict(excludeS[i]=>subs))
    end
    return excludeS
end