"""
Simple unimodal square function
True argmax - 0.25
"""
tf1(x:: Float64) = -(x - 0.25)^2 + 1 


"""
Sin wave with upward trend
True argmax - 0.787398
"""
tf2(x:: Float64) = sin(10 * x) + 0.2 * x


"""
Narrow guassian peak
True argmax - 0.7
"""
tf3(x:: Float64) = exp(-10000 * (x - 0.7)^2)


"""
2 Guassian peaks at different heights
True argmax - 0.487902
"""
tf4(x:: Float64) = exp(-50 * (x - 0.5)^2) + 0.8 * exp(-10 * (x - 0.075)^2)


"""
Kinked abs function with negative peak (non-differentiable)
True argmax - 0.7
"""
tf5(x:: Float64) = -abs(x - 0.7) - 1


"""
Linear function with max at corner
True argmax = 1.0
"""
tf6(x:: Float64) = x


"""
Guassian peak with sin bumps
True argmax - 0.501267
"""
tf7(x:: Float64) = exp(-10000 * (x - 0.5)^2) + 0.5 * sin(50 * x)


"""
Z modified so differnece in peaks is even smaller - can vary the denominator of x/10000 to see how it impacts grid search
True argmax = 0.9000005
"""
tf8(x:: Float64) = exp(-(10 * x - 1)^2) + exp(-(10 * x - 5)^2) + exp(-(10 * x - 9)^2) + x/10000


"""
Sin wave with guassian noise in middle and upward trend
True argmax = 0.518362
"""
tf9(x:: Float64) = sin(100x) * exp(-0.5(x - 0.5)^2) + 0.01x

"""
Linear increase that has a jump discontinuit followed by a local quadratic max
True argmax = 0.4
"""
function tf10(x::Float64) 
   if x < 0.4
        return x
    else
        return -(x-0.6)^2 - 0.2
    end
end

"""
Flat section that jumps to a quadratic curve
True argmax = 0.8
"""
function tf11(x::Float64)
    if x < 0.2
        return x
    elseif x <= 0.4
        return 0.2
    else
        return 1.0 - (x - 0.8)^2
    end
end

"""
Short strangle payoff but slightly domed
True argmax = 0.5 
"""
function tf12(x::Float64)
    if x < 0.3
        return x
    elseif x <= 0.7
        return 0.5 - 0.02*(x - 0.5)^2
    else
        return 1.0 - x
    end
end
"""
Function generator for testing robustness to small peaks
Denominator is adjustabel parameter -> as increases max peak should get harder to detect
Bounds: 0 < x1 < 0.5 & 0.5 < x2 < 1
True argmax = x2 (default pi/5)
"""
function make_test_func(denom:: Int; x1=pi/10, x2=pi/5, k=5000.0)
    function test_func(x:: Float64)
        peak1 = exp(-k * (x - x1)^2)
        peak2 = (1 + 1/denom) * exp(-k * (x - x2)^2)
        return peak1 + peak2
    end
    return test_func
end

"""
Function to generate 2 random points for x1 and x2 in the test_func. Allows robustness of placement of peaks
"""
function gen_random_points(; min_distance:: Float64 = 0.1)

    while true

        x1, x2 = sort(rand(2)) # sort ensures x1 is always smaller

        if x2 - x1 >= min_distance # implement the min distance to prevent peak inteference
            return x1, x2 # x1 < x2
        end
    end
end