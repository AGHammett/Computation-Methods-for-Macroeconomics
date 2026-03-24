using Plots
using Optim

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
    sigma = 2.0,
    T = 24.0,
    x = 8.0,
    epsilon = 1.0,
    alpha = 0.2,
    h_grid = range(0, 12, length = 500),
    h_eval = 10.0,
    y_lims = (16, 20.5)
)
    # initalise the basis of the plot. 
    p = plot(
        xlabel = "Hours worked (h)",
        ylabel = "Utility",
        title = "Z(h)",
        legend = :outerright,
        size = (1000, 600),
        linewidth = 2,
        ylims = y_lims,
        bottom_margin = 8Plots.mm,
        top_margin = 8Plots.mm,
        left_margin = 8Plots.mm,
        right_margin = 12Plots.mm,
        grid = true,
        gridalpha = 0.2
    )

    # find which value has been passed in to vary and loop through the arguments
    for v in values
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

        # evaluate Z at chose arguments over the grid
        Z_vals = [Z(w_i, h, a_i, g_i, s_i, T_i, x_i, e_i, A_i) for h in h_grid]
        plot!(p, h_grid, Z_vals, label = "$(param) = $(round(v, digits=2))", xlims = (0, 12)) # plot the function over the grid

        # plot the utility at h_eval for each value of the parameter passed in 
        if 0.0 <= h_eval <= T_i
            z_eval = Z(w_i, h_eval, a_i, g_i, s_i, T_i, x_i, e_i, A_i)
            plot!(p, [0.0, h_eval], [z_eval, z_eval],
                label = false,
                color = :black,
                linestyle = :dash,
                linewidth = 0.7,
                alpha = 0.8
            )
        end

        # find details for optimum marker
        res = optimize(h -> -Z(w_i, h, a_i, g_i, s_i, T_i, x_i, e_i, A_i), 0.0, T_i)
        h_star = Optim.minimizer(res)
        z_star = -Optim.minimum(res)

        # plot optimium
        if first(h_grid) <= h_star <= last(h_grid)
            scatter!(p, [h_star], [z_star],
                label = "optimisation solution",
                markershape = :xcross,
                markersize = 6,
                markerstrokewidth = 1
            )
        end
    end

    # plot the threshold value of x passed in
    vline!(p, [x],
        label = "Threshold x",
        color = :coral2,
        linestyle = :solid,
        linewidth = 1,
        alpha = 0.8
    )
    # plot the vertical line showing h_eval
    vline!(p, [h_eval],
        label = "Overtime ĥ = $(Int(round(h_eval)))",
        color = :orangered,
        linestyle = :solid,
        linewidth = 1
    )

    display(p)
    # save for use in document
    savefig(p, "q2plot.pdf")
    return p
end