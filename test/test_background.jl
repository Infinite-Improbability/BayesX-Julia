using Statistics

function test_uniform_background()
    @testset "Uniform Background" begin
        @test_throws AssertionError Uniform_Background(-1, 100)
    
        @testset "Checking length of spectrum" begin
            for len = [1, 10, 100]
                @test length(Uniform_Background(0, len)) == len
            end
        end

        @testset "Checking average rate correct" begin
            for rate = 1:1000:100
                spectrum = Uniform_Background(rate, 10000)
                avg = mean(spectrum)
                err = std(spectrum)
                @test avg - err <= rate <= avg + err
            end
        end    
    end
end

test_uniform_background()