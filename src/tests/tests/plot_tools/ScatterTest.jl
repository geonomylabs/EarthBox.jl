module ScatterTest

import EarthBox.PlotToolsManager: Scatter
import CairoMakie

function run_test()::Nothing
    println("Running ScatterTest")
    x_data = rand(100) * 10
    y_data = rand(100) * 10
    color_data = rand(100)

    fig = CairoMakie.Figure(size=(600, 400))
    ax = CairoMakie.Axis(fig[1, 1], title="Test Scatter Plot")

    # Test the plot_scatter function
    Scatter.plot_scatter(
        ax, 
        :viridis,  # colormap
        x_data, 
        y_data, 
        color_data,
        20.0,      # marker size
        0.0,       # min value
        1.0,       # max value
        [0.0, 0.25, 0.5, 0.75, 1.0],  # ticks
        label="Test Values"
    )

    # Save the test plot
    CairoMakie.save("test_scatter_plot.png", fig)
    println("Test scatter plot saved as 'test_scatter_plot.png'")

end

end # module