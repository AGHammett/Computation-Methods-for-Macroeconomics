using Plots
using Optim


# =============================== Functions that answer questions directly ==============================================

function q3_a()
    
    println("Running three peak implementation...")
    find_three_peaks(Z)
    println("-----------------------------------------------------")

    println("Running peak search implementation...")
    results = find_peaks(Z)
    z_star_max = -Inf
    h_star_max = 0
    for (i, result) in enumerate(results)
        h_star = Optim.minimizer(result)
        z_star = -Optim.minimum(result)
        
        println("Peak $i")
        println("H*: $h_star")
        println("Z*: $z_star")
        if z_star > z_star_max
            h_star_max = h_star
            z_star_max = z_star
        end
    end
    println("Max H*: $h_star_max")
    println("Max Z*: $z_star_max")

end

function q3_b(M:: Int = 100)

    println("Part 1 - evalutating h at m = 2 & m = 3:")
    h_2 = h(Z, 2)
    h_3 = h(Z, 3)
    println("h(2): $h_2")
    println("h(3): $h_3")


    println("Part 2 - graphing h over m:")
    m_range = 2 : 1 : M #set up array of xm starting from 2
     h_m_array = [h(Z, m) for m in m_range] # f can be any function that satisfies h

    p = plot(m_range, h_m_array, label = "h(m)") 
    hline!(p, [0.90005], label = "True Global Solution")
    display(p)
    savefig(p, "m1000.pdf")
    return p

end

function q3_c(f:: Function; grid_points = 1000)

    println("Optimising input function using Grid Search...")
    maximiser = global_solution(GridSearch(), f ; grid_points = grid_points)
    println("Global solution found by grid search: $maximiser")

    println("-----------------------------------------------------")

    println("Optimising input function using Peak Search...")
    maximiser = global_solution(PeakSearch(), f ; grid_points = grid_points)
    println("Global solution found by peak search: $maximiser")

end

# =============================== Functions that are used above ==============================================


"""
Z function as defined by assignment. Used to pass in as an argument to optimsiation or graphing functions
functional approach will allow for other functions with the same signature to be passed into the functions below
"""
function Z(h :: Float64):: Float64

    return exp(-(10 * h - 1)^2) + exp(-(10 * h - 5)^2) + exp(-(10 * h - 9)^2) + h/100

end

"""
Function taking one or multiple functions with a single argument and plots them on the same plot
"""
function plot_function(fs:: Function... ; grid_points:: Int = 100)

    grid = range(0, 1, length = grid_points)

    p = plot()
    for (i, f) in enumerate(fs) # loop through passed in functions and plot on the same object
        plot!(p, grid, f.(grid), label="Func $i")
    end

    display(p)
    return p
end

"""
Same function as above but saves plot with specified file name.
"""
function plot_function(export_name::String, fs:: Function...; grid_points:: Int = 100)

    grid = range(0, 1, length = grid_points)

    p = plot()
    for (i, f) in enumerate(fs) # loop through passed in functions and plot on the same object
        plot!(p, grid, f.(grid), label="Func $i")
    end

    display(p)
    savefig(p, export_name)
    return p
end

"""
Basic optimisation function to optimise a one dimensional function.
Optim's optimize function uses Brent by defualt when given a range
    
"""
function optimise_basic(f:: Function, lower_bound:: Float64, upper_bound:: Float64):: Tuple{Any, Float64, Float64}
    
    # optim optimisation can only minimise so we use the negative of the function to obtain a maximum
    result = optimize(x -> -f(x), lower_bound, upper_bound)  #Note optim usese a Z not S!!!
    h_star = Optim.minimizer(result) # optim.minizer will recover the argument that minises the function
    y_star = -Optim.minimum(result) # since we use the negative of the function we have to negative y_star
    
    return result, h_star, y_star

end

"""
Function to find the 3 local maxima of Z and the global max.
Elementary method of taking three subsets of the domain and optimising seperatly - by observing the graph it's obvious which value will work
"""
function find_three_peaks(f:: Function)

    # elementary function to find the three points based off a look at the graph. Gives correct result but isn't elegant at all

    ranges = ((0.0, 0.3), (0.3, 0.7), (0.7, 1.0)) # 1 max contained in each range based off viewing graph

    # create data containers
    results = []
    h_stars = Float64[]
    z_stars = Float64[]

    for range in ranges

        result, h_star, z_star = optimise_basic(f, range[1], range[2]) # apply the optimise_basic function from above

        # print and store results
        println("Range: $range")
        println("H*: $h_star")
        println("Z*: $z_star")

        push!(results, result)
        push!(h_stars, h_star)
        push!(z_stars, z_star)
    end

    global_max_index = argmax(z_stars) # to find the global maximum we look at the y values -> argmax lets us take its location in the array
    h_star_max = h_stars[global_max_index] # collect the h* and y* using the index
    z_star_max = z_stars[global_max_index]

    println("Max H*: $h_star_max")
    println("Max Y*: $z_star_max")

end

"""
Generates a grid of length h. Evaluates the function passed it at each grid points and returns the point that generated the max
"""
function h(f:: Function, m:: Int):: Float64

    #create an array with grid points
    grid = range(0.0, 1.0 , length =  m)
    outputs = f.(grid) # make use of Julia's function broadcasting to avoid looping.
 
    return grid[argmax(outputs)] # return the point that have the max value in outputs

end


# here I'm creating a new data type so use in the dispatch of the global_solution function to determine which algorithm it uses 
abstract type GlobalSolutionAlgorithm end # abstract type determines the class of the type but can't be instantiated
struct GridSearch <: GlobalSolutionAlgorithm end # <: inherits the type - no data assigned here so no memory needed
struct PeakSearch <: GlobalSolutionAlgorithm end


"""
Global optimiser with GridSearch
Uses a grid search follow by an Optim search method

f - one dimensional and over the domain [0, 1]
grid_point - int determining density of grid. Higher denisty -> greater accuracy but more expensive
search_method - can take either Brent() or GoldenSection() 

returns the maximiser
"""
function global_solution(::GridSearch, f::Function; grid_points::Int = 1000, search_method = Brent())::Float64
    grid = range(0.0, 1.0, length=grid_points)
    
    max_val = f(grid[1]) # evaluate first point to setup max_val 
    max_index = 1 # first point
    
    for i in 2:grid_points # in Julia looping actually gives better performance that f.(grid) due to memory efficiency
        val = f(grid[i])
        if val > max_val # only change point if it is highest seen
            max_val = val 
            max_index = i
        end
    end
    #grab indices either side to form the search interval 
    left_index = max(1, max_index - 1) # makes sure not to index out of bounds
    right_index = min(grid_points, max_index + 1)

    # optimise using the negative of the function - Optim only minimises so we transform 
    result = optimize(x -> -f(x), grid[left_index], grid[right_index], search_method)
    maximiser = Optim.minimizer(result)

    # check if there's a corner - brent/golden won't find it exactly
    for i in (0.0, 1.0)
        if f(i) > f(maximiser)
            maximiser = i
        end
    end
    return maximiser
end

#local points methods

"""
Function to find all maximums of a function over the domain [0,1] as long as they don't plateau

returns an array of Optim.optimize structs
"""
function find_peaks(f::Function; grid_points::Int = 1000, search_method = Brent())
    grid = range(0.0, 1.0, length=grid_points)
    peak_indices = Int[] # generate container for maxima
    
    #set up inital points for sliding window comparison
    y_left = f(grid[1])
    y_mid = f(grid[2])

    for i in 3:grid_points # loop through rest of array to compare all points
        y_right = f(grid[i])

        if y_mid > y_left && y_mid > y_right # points which are greater than the adjacent ones
            push!(peak_indices, i - 1) # add the index of the middle not the point
        end

        #update y values for next iteration
        y_left = y_mid
        y_mid = y_right
    end
    
    # if no peaks we may have an stricly increasing/decreasing function so we optimise around the whole interval
    if isempty(peak_indices)
        results = [optimize(x -> -f(x), 0.0, 1.0, search_method)] # needs to be a vector for consistent outputs
    else # otherwise use otim to find all peaks around the indices it's given
        results = [optimize(x -> -f(x), grid[i - 1], grid[i + 1], search_method) for i in peak_indices]
    end

    return results
end

"""
Global optimiser with PeakSearch
Finds all peaks in the inerval [0,1] and chooses the highest. It checks the highest peak against the edges
in case there was an interior peak but corner solution.

f - one dimensional and over the domain [0, 1]
grid_point - int determining density of grid. Higher denisty -> greater accuracy but more expensive
search_method - can take either Brent() or GoldenSection() 

returns the maximiser
"""
function global_solution(::PeakSearch, f:: Function ; grid_points:: Int = 1000,  search_method = Brent())::Float64
    results = find_peaks(f; grid_points = grid_points, search_method = search_method)

    results_minimums = Optim.minimum.(results) # use minimum here since results are from -f
    maximiser = Optim.minimizer(results[argmin(results_minimums)])
    # now must check against the ends for corner solutions
    for i in (0.0, 1.0)
        if f(i) > f(maximiser)
            maximiser = i
        end
    end
    return maximiser
end

# default dispatch of global_solution when tag is excluded
function global_solution(f:: Function ; grid_points:: Int = 1000,  search_method = Brent())
    return global_solution(GridSearch(), f ; grid_points = grid_points, search_method = search_method)
end

"""
Old version of the grid search optimiser
I thought vectorising the function would be slower but in Julia the memory overhead
    actually makes it slower than looping so this function is simply a worse version
"""
function global_solution_deprecated(f:: Function ; grid_points:: Int = 1000,  search_method = Brent()) # note kwargs allow flexibility

    # step 1 - grid search to identift potential maxima -> h can be reused here
    grid = (range(0.0, 1.0 , length =  grid_points))
    # step 2 - broadcast the function over the grid to get output at each point
    outputs = f.(grid)
    #step 3 - grab position of max point
    max_index = argmax(outputs)

    # create indices to get the bracket for testing - use min and max to cover corner cases.
    left_index = max(1, max_index - 1)
    right_index = min(grid_points, max_index + 1)

    # step 4 - optimise using the negative of the function - Optim only minimises so we transform 
    result = optimize(x -> -f(x), grid[left_index], grid[right_index], search_method)

    return result

end