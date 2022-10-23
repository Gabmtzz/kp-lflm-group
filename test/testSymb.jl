using  Symbolics

@testset "test Symbolic functions" begin
    @testset "test symbolic expression generation" begin
        @testset "test  simple aritmetic operations" begin
        
            basicOp =  ["+"; "-"; "*"; "/"; "^"]

            for operation  in basicOp
                
                strop= "1"*operation*"2"
                symbTestS = eval(Meta.parse(strop))
                SymbExprS = KPpack.StrtoSymbConv(strop)
                @test SymbExprS == symbTestS
            end

            @testset "test a complex aritmetic operation" begin
                exTest = (5-45+6+2)^6+(3/4+6)/12+34
                exQ = KPpack.StrtoSymbConv("(5-45+6+2)^6+(3/4+6)/12+34")
                @test exTest == exQ
             end
            
        end

        

        @testset "test other operations" begin

            othOp = ["sin"; "cos"; "tan"; "atan"; "sqrt"; "exp"]

            for operation in othOp
                strop = operation*"(1.424)"
                symbTestS = eval(Meta.parse(strop))
                SymbExprS = KPpack.StrtoSymbConv(strop)
                @test SymbExprS == symbTestS
            end
            
            @testset "test complex operation" begin
                exTest = sqrt(3+2)- exp(0.02)-(sin(2.2)^2)
                exQ = KPpack.StrtoSymbConv("sqrt(3+2)- exp(0.02)-(sin(2.2)^2)")
                @test exTest == exQ
            end
        end

        @testset "test expression with variables" begin
            a = Symbol("a"); a = @variables $a; a = a[1]
            b = Symbol("b"); b = @variables $b; b = b[1]
            c = Symbol("c"); c = @variables $c; c = c[1]

            testExpr = sin(a+b)+exp(a^3)/(3+sqrt(b))

            exprStrS = KPpack.StrtoSymbConv("sin(a+b)+exp(a^3)/(3+sqrt(b))")
            @test exprStrS - testExpr == 0

        end
    end
end

