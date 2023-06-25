#ifndef _elaspect_global_h
#define _elaspect_global_h

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace elaspect
{
  namespace constants
  {
    /**
     * Number of seconds in a year [s]
     */
    constexpr double year_in_seconds = 60*60*24*365.2425;

    /**
     * Zero degrees Celsius to Kelvin [K]
     */
    constexpr double celsius_to_kelvin = 273.15;

    /**
     * Gas constant (also known as R) [J K^-1 mol^-1]
     */
    constexpr double gas_constant = 8.3144621;
    /**
     * Avogadro's constant [mol^-1]
     */
    constexpr double avogadro = 6.02214129e23;
    /**
     * Gravitational constant [m^3 kg^-1 s^-2]
     */
    constexpr double big_g = 6.67384e-11;
  }

  using constants::year_in_seconds;

  /**
   * A typedef that denotes the BOOST stream type for reading data during
   * serialization. The type chosen here is a binary archive which we
   * subsequently will have to un-compress.
   */
  using iarchive = boost::archive::binary_iarchive;

  /**
   * A typedef that denotes the BOOST stream type for writing data during
   * serialization. The type chosen here is a binary archive which we compress
   * before writing it into a file.
   */
  using oarchive = boost::archive::binary_oarchive;

  class QuietException {};
}

/**
 * Print a header into the given stream that will be written both to screen
 * and to the log file and that provides basic information about what is
 * running and with how many processes.
 */
template <class Stream>
void print_elaspect_header(Stream &stream);

/**
 * A macro that is used in instantiating the elASPECT classes and functions for
 * both 2d and 3d. Call this macro with the name of another macro that when
 * called with a single integer argument instantiates the respective classes
 * in the given space dimension.
 */
#define ELASPECT_INSTANTIATE(INSTANTIATIONS) \
  INSTANTIATIONS(2) \
  INSTANTIATIONS(3)

#endif
