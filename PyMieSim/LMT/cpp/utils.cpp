double
GetSizeParameter(const double Diameter,
                 const double Wavelength,
                 const double nMedium)
{
  return PI * Diameter / (Wavelength / nMedium);
}
