include("q3_tfunctions.jl")
include("q3.jl")

using Optim
using BenchmarkTools

const TESTS = (
    (name = "tf1", f = tf1, target = 0.25),
    (name = "tf2", f = tf2, target = 0.787398),
    (name = "tf3", f = tf3, target = 0.7),
    (name = "tf4", f = tf4, target = 0.487902),
    (name = "tf5", f = tf5, target = 0.7),
    (name = "tf6", f = tf6, target = 1.0),
    (name = "tf7", f = tf7, target = 0.501267),
    (name = "tf8", f = tf8, target = 0.90005),
    (name = "tf9", f = tf9, target = 0.518362),
    (name = "Z",   f = Z,   target = 0.90005))

function plot_test_functions()

    plot_function(tf1, tf2, tf3, tf4, tf5, tf6, tf7, tf8, tf9, grid_points = 1000)

end

function test_optimisation(grid_points::Int = 1000)

    for (i, (name, f, target)) in enumerate(TESTS)

        println("Testing $name")
        result_b = global_solution(f, grid_points = grid_points)
        result_g = global_solution(f, grid_points = grid_points, search_method = GoldenSection())
        println("Golden Section Difference: ", abs(Optim.minimizer(result_b) - target)) 
        println("Golden Section Difference: ", abs(Optim.minimizer(result_g) - target)) 

    end

end

function test_optimisation_speed()

    grid_points_test = (10, 100, 1000, 10000, 100000, 1000000, 10000000)

    b_times = Float64[]
    g_times = Float64[]

    for grid_points in grid_points_test
    # use tf8 for the benchmark test - all so similar
        b_time = @benchmark global_solution($tf8, grid_points = $grid_points)
        g_time = @benchmark global_solution($tf8, grid_points = $grid_points, search_method = GoldenSection())

        push!(b_times, median(b_time).time / 1e6) 
        push!(g_times, median(g_time).time / 1e6)

    end

    println(b_times)

    p_log = plot(collect(grid_points_test), b_times, 
             label = "Brent", xscale = :log10, yscale = :log10,
             xlabel = "Grid Points", ylabel = "Time (ms)", 
             title = "Optimiser Speed vs Grid Size (Log Scale)",
             xticks = collect(grid_points_test),
             yticks = ([0.001, 0.01, 0.1, 1, 10, 100, 1000], ["0.001", "0.01", "0.1", "1", "10", "100", "1000"]))
    plot!(collect(grid_points_test), g_times, label = "Golden Section")

    p_lin = plot(collect(grid_points_test), b_times, 
             label = "Brent",
             xlabel = "Grid Points", ylabel = "Time (ms)", 
             title = "Optimiser Speed vs Grid Size (Linear Scale)",
             xticks = (collect(grid_points_test), ["10", "100", "1K", "10K", "100K", "1M", "10M"]),
             xrotation = 45)
    plot!(collect(grid_points_test), g_times, label = "Golden Section")
    
    display(p_log)
    display(p_lin)

    savefig(p_log, "optimiser_speed_log.png")
    savefig(p_lin, "optimiser_speed_linear.png")

    return p_log, p_lin

end

"""
Time of each of the test functions - all take 1.600ns so arbitrary which is used for testing
"""
function test_function_times()

    for (name, f, target) in TESTS
        function_time = @benchmark $f(0.5)
        println(name, median(function_time))
    end
end

function test_ngrid_accuracy(;tolerance = 1.0e-3)

    grid_points = 1000 # start at 1000 since I know this passes all tests
    failed_b = Dict{String, Int}()
    failed_g = Dict{String, Int}()
    
    while length(failed_b) < length(TESTS) && length(failed_g) < length(TESTS)
        
        for (name, f, target) in TESTS
            
            haskey(failed_b, name) && haskey(failed_g, name) && continue

            result_b = global_solution(f, grid_points = grid_points)
            result_g = global_solution(f, grid_points = grid_points, search_method = GoldenSection())

            if abs(Optim.minimizer(result_b) - target) > tolerance
                failed_b[name] = grid_points
            end

            if abs(Optim.minimizer(result_g) - target) > tolerance
                failed_g[name] = grid_points
            end

        end

        grid_points -= 1
        # need a break statement because some functions work with 2 points and 1 causes range to throw an error
        if grid_points == 1
            break
        end
    end

    return failed_b, failed_g

end

function test_close_peaks(;tolerance = 1.0e-3)

    denominator = 10
    grid_points = 100
    denoms = Int[]
    succ_grid_points = Int[]

    while denominator < 1.01e8

        println("Testing $denominator")

        f = make_test_func(denominator, k=2e5)

        while grid_points < 1e8

            result = global_solution(f, grid_points = grid_points)

            if abs(Optim.minimizer(result) - (pi / 5)) > tolerance
                grid_points *= 2
            else
                push!(denoms, denominator)
                push!(succ_grid_points, grid_points)
                break
            end
        end

        denominator *= 10
    end
    return denoms, succ_grid_points
end
