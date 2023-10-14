"""
Scattering efficiency of a sphere
=================================

PyMieSim makes it easy to create a source and a scatterer. With these objects
defined, it is possible to use PyMieSim to find the scattering efficiency of the
scatterer. This feature can be used to plot a graph of the scattering efficiency
of a sphere as a function of the permittivity and the size parameter.
"""

# sphinx_gallery_thumbnail_path = '../images/Extras/ScatteringEfficiency.png'



def run():
    from PyMieSim.source import PlaneWave
    from PyMieSim.scatterer import Sphere

    # import other librairies
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib

    # set the source using PyMieSim
    Source = PlaneWave(
        wavelength=450e-9,
        polarization=0,
        amplitude=1
    )

    # create an empty list that will contain the heatmap values
    heatmap = []

    # loop through all the values to generate the heatmap
    # first, loop through the diameter values
    for i in range(1, 517):
        x = []
        # second, loop through the index values
        for j in range(-100, 501):
            if j == 0:
                continue

            Scat = Sphere(
                diameter=i / 3 * 1e-9,
                source=Source,
                index=(j / 10)**0.5
            )  # square root since permittivity, is the index squared

            # Get the scattering efficiency using GetProperties()
            prop = Scat.GetProperties()
            Qsca = prop[1]
            x.append(Qsca)
        # append list (row) to the heatmap (list of lists)
        heatmap.append(x)

    # convert the heatmap to a numpy array
    data = np.array(heatmap)

    # create the plot
    fig, ax = plt.subplots()
    # show the data and adjust the color scale
    im = ax.imshow(data, norm=matplotlib.colors.LogNorm(vmin=0.1, vmax=10), cmap='viridis')

    # graph title
    ax.set_title("Scattering efficiency of a sphere")

    # x axis settings
    ax.set_xlabel("Permittivity")
    ax.set_xticks(np.linspace(0, len(heatmap[0]), 7))
    ax.set_xticklabels([-10, 0, 10, 20, 30, 40, 50])

    # y axis settings
    ax.set_ylabel("Size parameter")
    ax.invert_yaxis()
    ax.set_yticks(np.linspace(0, len(heatmap), 7))
    ax.set_yticklabels([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])

    # colorbar settings
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)

    # display the plot in a tight layout
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    run()

# -

# -
