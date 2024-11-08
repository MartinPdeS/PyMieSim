import io
from flask import Flask, render_template, request, send_file
import matplotlib.pyplot as plt
import numpy as np

app = Flask(__name__)

def plot_mie_scattering(diameter, wavelength, refractive_index):
    """
    Generate a plot based on diameter, wavelength, and refractive index.

    Parameters
    ----------
    diameter : float
        Diameter of the particle.
    wavelength : float
        Wavelength of the incident light.
    refractive_index : float
        Refractive index of the particle.
    """
    # Generate x values (for simplicity, we'll use angles)
    x = np.linspace(0, np.pi, 100)

    # Simulate a simple Mie scattering pattern (for demonstration)
    intensity = np.abs(np.sin(diameter * x / wavelength) * np.exp(-refractive_index * x))

    # Plot
    plt.figure(figsize=(8, 4))
    plt.plot(x, intensity, label=f'Diameter={diameter}µm, Wavelength={wavelength}nm, n={refractive_index}')
    plt.title("Simulated Mie Scattering Intensity")
    plt.xlabel("Angle (radians)")
    plt.ylabel("Intensity")
    plt.legend()
    plt.grid(True)

    # Save the plot to a BytesIO object
    img = io.BytesIO()
    plt.savefig(img, format='png')
    img.seek(0)
    plt.close()  # Close the plot to free memory
    return img

@app.route('/', methods=['GET', 'POST'])
def index():
    """
    Display the form and show the plot.
    """
    if request.method == 'POST':
        # Get parameters from the form
        diameter = float(request.form['diameter'])
        wavelength = float(request.form['wavelength'])
        refractive_index = float(request.form['refractive_index'])

        # Generate the plot image
        img = plot_mie_scattering(diameter, wavelength, refractive_index)
        return send_file(img, mimetype='image/png')
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
