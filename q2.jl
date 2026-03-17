using Plots

U(c:: Real, gamma:: Real) = (c ^ (1 - gamma)) / (1 - gamma)

V(h:: Real, sigma:: Real, T::Real) = ((T - h) ^ (1 - 1 / sigma)) / (1 - 1 / sigma)

W(h:: Real, x:: Real, epsilon:: Real, alpha:: Real) = alpha * ((1 + max(0, h - x)) ^ epsilon - 1)

con(w:: Real, h::Real, a::Real) = w * h + a 

Z(w, h, a, gamma, sigma, T, x, epsilon, alpha) = U(con(w, h, a), gamma) + V(h, sigma, T) - W(h, x, epsilon, alpha)

function plot_utility_variation(; 
    param::Symbol,
    values,
    w = 10.0,
    a = 0.5,
    gamma = 0.8,
    sigma = 2,
    T = 24.0,
    x = 8.0,
    epsilon = 1.0,
    alpha = 0.2,
    h_grid = range(0, 12, length = 500)
)

    p = plot(
        xlabel = "Hours worked (h)",
        ylabel = "Utility",
        title = "Utility as function of hours worked")
    
    Z_base = [Z(w, h, a, gamma, sigma, T, x, epsilon, 0.0) for h in h_grid]
    #plot!(p, h_grid, Z_base, label = "alpha = 0.0", linestyle = :dash, linewidth = 1, linecolor = "blue")

    for v in values

        # baseline parameters
        w_i, a_i, g_i, s_i, T_i, x_i, e_i, A_i = w, a, gamma, sigma, T, x, epsilon, alpha

        if param == :w
            w_i = v
        elseif param == :a
            a_i = v
        elseif param == :gamma
            g_i = v
        elseif param == :sigma
            s_i = v
        elseif param == :T
            T_i = v
        elseif param == :x
            x_i = v
        elseif param == :epsilon
            e_i = v
        elseif param == :alpha
            A_i = v
        else
            error("Unknown parameter")
        end

        Z_vals = [Z(w_i, h, a_i, g_i, s_i, T_i, x_i, e_i, A_i) for h in h_grid]

        plot!(p, h_grid, Z_vals, label = "$(param) = $(v)")
    end

    display(p)
end