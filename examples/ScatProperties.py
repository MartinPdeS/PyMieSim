def run():
    from PyMieSim.Scatterer import Sphere
    from PyMieSim.Source    import PlaneWave

    Source = PlaneWave(Wavelength   = 450e-9,
                       Polarization = 0,
                       E0           = 1)

    Scat = Sphere(Diameter    = 800e-9,
                  Source      = Source,
                  Index       = 1.4)

    print(Scat.Properties)


if __name__ == '__main__':
    run()


# _____________________OUTPUT_____________________
#       Object:          Dictionary
#        Keys:            Efficiencies, cross-sections, others
#        Structured data: Yes
#        Method:          <Plot>
#        Shape:           [7, 1]
#         ========================================
# --------------------------------------------------
# Efficiencies   | Qsca        | 4.029799032677242
# --------------------------------------------------
# Efficiencies   | Qext        | 4.029799032677242
# --------------------------------------------------
# Efficiencies   | Qabs        | 0.0
# --------------------------------------------------
# Efficiencies   | Qback       | 4.973830378796597
# --------------------------------------------------
# Efficiencies   | Qratio      | 1.2342626365395144
# --------------------------------------------------
# Efficiencies   | Qpr         | 0.7925897259835781
# --------------------------------------------------
# cross-sections | Csca        | 2.03e+00 μm²
# --------------------------------------------------
# cross-sections | Cext        | 2.03e+00 μm²
# --------------------------------------------------
# cross-sections | Cabs        | 0 m²
# --------------------------------------------------
# cross-sections | Cback       | 2.50e+00 μm²
# --------------------------------------------------
# cross-sections | Cratio      | 6.20e+05 nano-m²
# --------------------------------------------------
# cross-sections | Cpr         | 3.98e+05 nano-m²
# --------------------------------------------------
# others         | area        | 5.03e+05 nano-m²
# --------------------------------------------------
# others         | index       | 1.4
# --------------------------------------------------
# others         | g           | 0.7925897259835781
# --------------------------------------------------
