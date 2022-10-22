@testset "test reading Params " begin
    arr=["GaAs" 6.98;
        "AlAs" 3.76;
        "InAs" 19.7;
        "GaSb" 13.4; 
        "HgTe" 4.1;
        "CdTe" 1.47]
    
    @testset "Test bulk Material properties for: $(arr[i,1])" for i ∈ 1:size(arr)[1]
        T = 200

        testparM = KPpack.Materials(arr[i,1],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        KPpack.params(testparM)
        α,β,alat = KPpack.TempPar(testparM.material)
        EgRef,alRef = (testparM.Eg - (α*T^2)/(T+β)), (testparM.al+alat*(T-300))
        EgTtest, alTtest = KPpack.aEgTemp(T,testparM.material,testparM.Eg,testparM.al)

        @test testparM.g1 == arr[i,2]
        @test EgTtest == EgRef
        @test alTtest == alRef
    end

    testparM = KPpack.Materials(arr[1,1],0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    KPpack.params(testparM)
    @testset "Test calculus F" begin
        
        S = (1/testparM.me)-(testparM.Ep*(testparM.Eg+(2/3)*testparM.delta)/(testparM.Eg*(testparM.Eg+testparM.delta)))
        Ftest = (S-1)/2

        KPpack.getF(testparM,"calc")
        @test testparM.F == Ftest
    end
    
    @testset "test Alloy properties" begin
        @testset "test the name creation" begin
            matAlloy = "AlGaAs_0.2"
            mat1Test = "AlAs"; mat2Test="GaAs"; compTest=0.2

            mat1,mat2,comp= KPpack.DetMat(matAlloy)
            @test mat1==mat1Test
            @test mat2==mat2Test
            @test comp==compTest
        end

        @testset "test alloy parameters" begin
            matAlloy = "AlGaAs_0.2"
            mat1Test = "AlAs"; mat2Test="GaAs"; compTest=0.2

            T = 0; opt = "calc"
            @testset "test function ParMat for bulk" begin
                Egtest = 1.519
                bulkStrct  = KPpack.Materials(mat2Test,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
                KPpack.ParMat(bulkStrct,T,opt)
                @test bulkStrct.Eg == Egtest
            end
            @testset "test function ParMat for alloy" begin
                Egtest = 1.7954
                alloyStrct  = KPpack.Materials(matAlloy,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
                KPpack.ParMat(alloyStrct,T,opt)
                @test alloyStrct.Eg ≈  Egtest
            end
        end
    end
end