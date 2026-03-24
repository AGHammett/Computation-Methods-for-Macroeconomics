include("q3_tfunctions.jl")
include("q3.jl")

using Optim
using BenchmarkTools
using DataFrames
using CSV

const TESTS = ( # test functions form q3_testfunctions
    (name = "tf1", f = tf1, target = 0.25),
    (name = "tf2", f = tf2, target = 0.787398),
    (name = "tf3", f = tf3, target = 0.7),
    (name = "tf4", f = tf4, target = 0.487902),
    (name = "tf5", f = tf5, target = 0.7),
    (name = "tf6", f = tf6, target = 1.0),
    (name = "tf7", f = tf7, target = 0.501267),
    (name = "tf8", f = tf8, target = 0.9000005),
    (name = "tf9", f = tf9, target = 0.518362),
    (name = "tf10", f = tf10, target = 0.4),
    (name = "tf11", f = tf11, target = 0.8),
    (name = "tf12", f = tf12, target = 0.5),
    (name = "Z",   f = Z,   target = 0.90005))

"""
Test optimiser against the tuple of test functions using Brent and GoldenSection
If verbose will print the difference for every test 
If not verbose will only notify on failure
"""
function test_optimisation(grid_points::Int = 1000; algorithm::GlobalSolutionAlgorithm = GridSearch(), verbose::Bool = false, tolerance = 1e-6)

    for (i, (name, f, target)) in enumerate(TESTS)

        println("Testing $name")
        # get optim results for each method
        result_b = global_solution(algorithm, f, grid_points = grid_points) # recall global_solution already takes negative of function so we pass in just the function
        result_g = global_solution(algorithm, f, grid_points = grid_points, search_method = GoldenSection())

        # compare maximiser to target - difference should be 0 for success
        b_diff = abs(result_b - target)
        g_diff = abs(result_g - target)
        
        if b_diff > tolerance || g_diff > tolerance # test function fails if outside of tolerance.
            println("Failed function ", name)
        end
        
        if verbose == true # verbose function will display exact difference
            println("Brent Difference: ", b_diff) 
            println("Golden Section Difference: ", g_diff) 
        end
    end

end

"""
Benchmark time and memory performance for the alternate algroithms in global_solution.
grid_points argument controls grid density used by solution methods 
Returns a Dict with algorithm as key and tuple of results (time, memory, allocs)
"""
function benchmark_optimisation_methods(grid_points::Int)
    results = Dict()

    for (name, f, _) in TESTS
        method_results = Dict{Symbol, Tuple{Float64, Int, Int}}()

        t1 = @benchmark global_solution(GridSearch(), $f; grid_points = $grid_points)
        est1 = mean(t1)
        method_results[:GridSearch] = (est1.time, est1.memory, est1.allocs)

        t2 = @benchmark global_solution(PeakSearch(), $f; grid_points = $grid_points)
        est2 = mean(t2)
        method_results[:PeakSearch] = (est2.time, est2.memory, est2.allocs)

        results[name] = method_results
    end

    return results
end

"""
Outputs benchmark results for easy comparison between methods
"""
function save_benchmark_results()

    # function to take in a result dict from the benchmark_optimisation_methods function and create a dataframe of them
    function results_to_dataframe(results::Dict, grid_points::Int)
        rows = []

        for (fname, methods) in results
            for (method, (time, memory, allocs)) in methods
                push!(rows, (
                    function_name = fname,
                    method = String(method),
                    grid_points = grid_points,
                    time_ns = time,
                    memory_bytes = memory,
                    allocs = allocs
                ))
            end
        end
        return DataFrame(rows)
    end

    # output results of 1,000 and 100,000 grid point tests
    df1 = results_to_dataframe(benchmark_optimisation_methods(1000), 1000)
    df2 = results_to_dataframe(benchmark_optimisation_methods(100000), 100000)

    df = vcat(df1, df2) # concat dfs
    CSV.write("benchmark_results.csv", df)

end

"""
graphs how long it takes to run optimiser at different grid densities
Due to benchmarking this function takes a while to run
"""
function plot_optimisation_speed()

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

"""
function that starts at a grid of 1000 and decreases resolution to see when optimisation fails 
    to find the global maximum. 
Test using both Brent and GoldenSearch for completion but it is technically redunadnat as failure
    arises due to grid inaccuracy -> this shows both methods give the same result
"""
function test_ngrid_accuracy(;tolerance = 1.0e-3)

    grid_points = 1000 # start at 1000 since I know this passes all tests

    # setup dictionaries for failed tests
    failed_b = Dict{String, Int}()
    failed_g = Dict{String, Int}()
    
    while length(failed_b) < length(TESTS) && length(failed_g) < length(TESTS)
        
        for (name, f, target) in TESTS
            
            # while the test isn't part of the failed dictionary try it. if not it will skip
            if !haskey(failed_b, name)
                result_b = global_solution(f, grid_points = grid_points)
                if abs(result_b - target) > tolerance # only add once out of tolerance
                    failed_b[name] = grid_points
                end
            end
            #same logic but for goldensearch
            if !haskey(failed_g, name)
                result_g = global_solution(f, grid_points = grid_points, search_method = GoldenSection())
                if abs(result_g - target) > tolerance
                    failed_g[name] = grid_points
                end
            end

        end

        grid_points -= 1
        # need a break statement because some functions work with 2 points and 1 causes range to throw an error
        if grid_points == 1
            break
        end
    end
    
    # print resutls
    println("Brent Method:")
    
    if isempty(failed_b)
        println("  All functions passed down to grid_points = $grid_points")
    else
        for (name, gp) in failed_b
            println("  $name first failed at grid_points = $gp")
        end
    end

    println("\nGolden Section Method:")
    if isempty(failed_g)
        println("  All functions passed down to grid_points = $grid_points")
    else
        for (name, gp) in failed_g
            println("  $name first failed at grid_points = $gp")
        end
    end
end


"""
Generates a function with 2 guassian peaks - one slightly higher than the other. Decresases size difference to test for failure. If failure it increase the grid density.
The default version tests where maxes are at pi/10 and pi/5
The random version will generate 2 points randomly for the peaks so is more adversarial - it can't register coincidental points. The optimiser needs to pass all random variations to move
    to the next difficulty level
At the end it prints how many grid points are needed for each difficulty level

Difficulty and grid density scaling rates are controlable parameters 
"""
function test_close_peaks(;algorithm:: GlobalSolutionAlgorithm = GridSearch(), tolerance = 1.0e-3, difficulty_scaler:: Real = 10, grid_scaler:: Real = 2, k = 5000.0,  random_tests::Int = 0)

    denominator = 10
    grid_points = 10 # note grid points never resets after successful test - harder tests will required at least as many
    denoms = Int[]
    succ_grid_points = Int[]

    while denominator < 1.01e8
        println("Testing $denominator")

        success = false

        # Repeat until current denominator passes
        while !success && grid_points < 1e8

            # randomise the points being tested against - more robust
            if random_tests > 0
                random_fails = 0

                for i in 1:random_tests
                    x1, x2 = gen_random_points()
                    f = make_test_func(denominator, x1=x1, x2=x2, k = k) 
                    target = x2 # x2 will be the bump that has a small difference added to it 
                    result = global_solution(algorithm, f, grid_points=grid_points)

                    if abs(result - target) > tolerance
                        random_fails += 1 
                        # double grid points and restart the test for all sets of points
                        grid_points *= grid_scaler
                        break
                    end
                end

                success = (random_fails == 0) # success stat only becomes true when all tests are passed

            # use default peaks in test function
            else
                f = make_test_func(denominator, k = k)
                target = pi / 5 # default target
                result = global_solution(algorithm, f, grid_points=grid_points)

                if abs(result - target) <= tolerance
                    success = true # true if it passes the single test here
                else
                    grid_points *= grid_scaler
                end
            end
        end  # ends while !success when success become true to move onto next difficult

        # Record result
        push!(denoms, denominator)
        push!(succ_grid_points, grid_points)

        # Move to next difficulty
        denominator *= difficulty_scaler
    end
    return denoms, succ_grid_points
end
