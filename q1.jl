using Plots
using Optim
using Roots

u(c::Float64, gamma:: Real)::Float64 = (c^(1 - gamma)) / (1 - gamma)

v(h::Float64, sigma:: Real)::Float64 = ((1 - h)^(1 - 1 / sigma)) / (1 - 1 / sigma)

function Z1(h::Float64 ; w::Real = 2, gamma:: Real = 3, sigma:: Real = 0.5, a:: Real = 0.5):: Float64

    c = w * h + a

    return u(c, gamma) + v(h, sigma)

end

function plot_function(f:: Function; p_range:: Tuple{Real, Real} = (0.0, 0.9), grid_points:: Int = 100)

    grid = range(p_range[1], p_range[2], grid_points)
    
    p = plot(grid, f.(grid), label = "Function")
    display(p)
    return p

end

function optimise_Z(;w::Real = 2, gamma:: Real = 3, sigma:: Real = 0.5, a:: Real = 0.5)

    result = optimize(h -> -Z1(h; w = w, gamma = gamma, sigma = sigma, a = a), 0.0, 0.9)
    return Optim.minimizer(result)
end

function find_a_min(grid_end:: Real, grid_size:: Int; tolerance = 1e-4)

    a_grid = range(0, grid_end, grid_size)

    for a in a_grid
        h = optimise_Z(a = a)
        if h < tolerance
            a_estimate = a
            break
        end
    end

end