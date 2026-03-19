using Roots
using Plots
using Optim


# note optim's default method is Nelder-Mead

"""
Z function as defined by assignment. Used to pass in as an argument to optimsiation or graphing functions
functional approach will allow for other functions with the same signature to be passed into the functions below
"""
function Z(h :: Float64):: Float64

    return exp(-(10 * h - 1)^2) + exp(-(10 * h - 5)^2) + exp(-(10 * h - 9)^2) + h/100

end

"""
Plots any function that takes a single Float64 argument over the interval [0,1]
Displays and returns plot
"""
function plot_function(f:: Function; grid_points:: Int = 100)


    # default of 100 points shows most functions smoothly - increase for better resolution with more detailed functions
    grid = range(0, 1, length = grid_points) # note range returns a range type NOT an array -> more memory effecient for looping & grids

    p = plot(grid, f.(grid), label = "Func")
    display(p)
    return p

end

"""
Alternate function taking multiple fucntions as specified above and plots them on the same plot
Uses multiple dispatch
"""
function plot_function(fs:: Function... ; grid_points:: Int = 100)

    grid = range(0, 1, length = grid_points)

    p = plot()
    for (i, f) in enumerate(fs)
        plot!(p, grid, f.(grid), label="Func $i")
    end

    display(p)
    return p
end

"""
Basic optimisation function to optimise a one dimensional function.
Optim's optimize function uses Brent by defualt when given a range
    
"""
function optimise_basic(f:: Function, lower_bound:: Float64, upper_bound:: Float64):: Tuple{Any, Float64, Float64}
    

    # optim optimisation can only minimise so we use the negative of the function to obtain a maximum
    result = optimize(x -> -f(x), lower_bound, upper_bound)  #Note optim usese a Z not S!!!
    h_star = Optim.minimizer(result)
    y_star = -Optim.minimum(result) # since we use the negative of the function we have to negative y_star
    
    return result, h_star, y_star

end

function find_three_peaks(f:: Function)

    # elementary function to find the three points based off a look at the graph. Gives correct result but isn't elegant at all

    ranges = ((0.0, 0.3), (0.3, 0.7), (0.7, 1.0))

    results = []
    h_stars = Float64[]
    y_stars = Float64[]

    for range in ranges

        result, h_star, y_star = optimise_basic(f, range[1], range[2])

        println("Range: $range")
        println("H*: $h_star")
        println("Y*: $y_star")

        push!(results, result)
        push!(h_stars, h_star)
        push!(y_stars, y_star)
    end

    global_max_index = argmax(y_stars) # to find the global maximum we look at the y values -> argmax lets us take its location in the array
    h_star_max = h_stars[global_max_index]
    y_star_max = y_stars[global_max_index]

    println("Max H*: $h_star_max")
    println("Max Y*: $y_star_max")

end

function h(f:: Function, m:: Int):: Float64

    #create an array with grid points
    grid = range(0.0, 1.0 , length =  m)
    outputs = f.(grid) # make use of Julia's function broadcasting to avoid looping.

    max_index = argmax(outputs)
    return grid[max_index] # can ignore max_index entirely here 

end

function q3_b(f:: Function, M:: Int)

    m_range = 2 : 1 : M #set up array of xm starting from 2
     h_m_array = [h(f, m) for m in m_range] # f can be any function that satisfies h

    p = plot(m_range, h_m_array, label = "h(m)") 
    hline!(p, [0.90005])
    display(p)
    return p

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

returns an Optim.optimizer struct
"""
function global_solution(::GridSearch, f:: Function ; grid_points:: Int = 1000,  search_method = Brent()) # note kwargs allow flexibility

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

#local points methods

"""
Function to find all maximums of a function over the domain [0,1] as long as they don't plateau

    Returns an array of Optim.optimize structs
"""
function find_peaks(f::Function; grid_points::Int = 1000, search_method = Brent())
    grid = range(0.0, 1.0, length=grid_points)
    outputs = f.(grid)

    peak_indices = Int[] # generate container for maxima
    
    for i in 2:grid_points-1 # not loop starts and ends 1 away form the edges or we'd have index errors
        if outputs[i] > outputs[i - 1] && outputs[i] > outputs[i + 1] # points which are greater than the adjacent ones
            push!(peak_indices, i) # add the index not the point
        end
    end
    
    # if no peaks we have a unimodal function so it can be handled directly
    if isempty(peak_indices)
        results = [optimize(x -> -f(x), 0.0, 1.0, search_method)] # needs to be a vector for consistent outputs
    else # otherwise use otim to find all peaks around the indices it's given
        results = [optimize(x -> -f(x), grid[i - 1], grid[i + 1], search_method) for i in peak_indices]
    end

    return results
end

"""
Global optimiser with PeakSearch
Finds all peaks in the inerval [0,1] and chooses the highest

f - one dimensional and over the domain [0, 1]
grid_point - int determining density of grid. Higher denisty -> greater accuracy but more expensive
search_method - can take either Brent() or GoldenSection() 

returns n Optim.optimizer struct
"""
function global_solution(::PeakSearch, f:: Function ; grid_points:: Int = 1000,  search_method = Brent())
    results = find_peaks(f; grid_points = grid_points, search_method = search_method)

    results_minimums = Optim.minimum.(results) # use minimum here since results are from -f
    return results[argmin(results_minimums)]
end

# default dispatch of global_solution when tag is excluded
function global_solution(f:: Function ; grid_points:: Int = 1000,  search_method = Brent())
    return global_solution(GridSearch(), f ; grid_points = grid_points, search_method = search_method)
end