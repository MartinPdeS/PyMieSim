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

extern const std::unordered_map<std::string, std::array<BesselErrors, 6>> errorMessages;
