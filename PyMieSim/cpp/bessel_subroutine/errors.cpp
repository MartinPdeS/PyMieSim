#pragma once

#include <string>
#include <array>
#include <unordered_map>

// Error handling main code inspired by this:
// https://stackoverflow.com/questions/47841783/is-there-any-advantage-in-using-stdoptional-to-avoid-default-arguments-in-a-fu
enum ErrorCode
{
  Success,
  InputError,
  Overflow,
  PartialLossOfSignificance,
  FullLossOfSignificance,
  AlgorithmTermination
};

typedef struct BesselErrors
{
  ErrorCode   errorCode;    // Equivalent of IERR in the original FORTRAN.
  std::string errorMessage; // Error message associated with the error.
} BesselErrors;

static const std::unordered_map<std::string, std::array<BesselErrors, 6>> errorMessages = {
  { "besselJ",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, Im(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) or order large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) or order too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "besselY",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, order is too large or abs(z) is too small or both" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) or order large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) or order too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "besselI",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, Re(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) or order large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) or order too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "besselK",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, order is too large or abs(z) is too small or both" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) or order large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) or order too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "hankelH",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, order is too large or abs(z) too small or both" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) or order large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) or order too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "airy",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, Re(2/3*z*sqrt(z)) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } }
  },
  { "biry",
    std::array<BesselErrors, 6>{
      BesselErrors{ Success, "Normal return\t -- Computation completed." },
      BesselErrors{ InputError, "Input error\t -- No computation" },
      BesselErrors{ Overflow, "Overflow\t -- No computation, Re(z) too large for scale=false" },
      BesselErrors{ PartialLossOfSignificance, "Loss of significance\t -- abs(z) large \n computation done but losses of significance by argument reduction produce less than half of machine accuracy" },
      BesselErrors{ FullLossOfSignificance, "Loss of significance\t -- abs(z) too large\n no computation because of complete loss of significance by argument reduction" },
      BesselErrors{ AlgorithmTermination, "Error\t -- no computation, algorithm termination condition not met" } } }
};
