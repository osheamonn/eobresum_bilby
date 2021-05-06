"""

Contains wrapper around TEOBResum code which satisfies the interace requirements of Bilby.

Assumes that EOB resum has is installed in the current python environment, and libconfig must also be available.

This file doesn't do everything EOBResum can, for instance it's possible to return all of the modes along with the strain (using 'output_hpc' keyword), but we just return the full strain.

EOBResum takes several keyword arguments which are passed along:

'use_mode_lm': Array of which modes to use. The model uses a 1D representation, k, of (l,m) modes with m>0 in all cases such that increasing k corresponds to increasing l. The conversion formula is k = l*(l-1)/2 + m -2.

E.g. 'use_mode_lm': [0, 1, 2]
will use the (2,1), (2,2) and (3,1) modes.

'ode_abstol': Absolute error tolerance for ode integrator
'ode_reltol': Relative error tolerance for ode integrator

"""

import numpy as np
import EOBRun_module
from scipy.fft import fft
from pycbc.types import TimeSeries
from pycbc.waveform.utils import taper_timeseries


def eccentric_binary_black_hole_eob_resum_nonspinning(
    frequency_array,
    mass_1,
    mass_2,
    eccentricity,
    chi_1,
    chi_2,
    luminosity_distance,
    theta_jn,
    phase,
    **kwargs
):
    domain = None
    if "FrequencyDomain" in kwargs and not "TimeDomain" in kwargs:
        domain = 0
        return native_frequency_domain(
            frequency_array,
            mass_1,
            mass_2,
            eccentricity,
            0.,
            0.,
            luminosity_distance,
            theta_jn,
            phase,
            **kwargs
        )
    elif "TimeDomain" in kwargs:
        domain = 1
        return fourier_transform_time_domain(
            frequency_array,
            mass_1,
            mass_2,
            eccentricity,
            0.,
            0.,
            luminosity_distance,
            theta_jn,
            phase,
            **kwargs
        )
    else:
        raise RuntimeError("Either TimeDomain or FrequencyDomain should be true")


def eccentric_binary_black_hole_eob_resum_aligned_spins(
    frequency_array,
    mass_1,
    mass_2,
    eccentricity,
    chi_1,
    chi_2,
    luminosity_distance,
    theta_jn,
    phase,
    **kwargs
):
    domain = None
    if "FrequencyDomain" in kwargs and not "TimeDomain" in kwargs:
        domain = 0
        return native_frequency_domain(
            frequency_array,
            mass_1,
            mass_2,
            eccentricity,
            chi_1,
            chi_2,
            luminosity_distance,
            theta_jn,
            phase,
            **kwargs
        )
    elif "TimeDomain" in kwargs:
        domain = 1
        return fourier_transform_time_domain(
            frequency_array,
            mass_1,
            mass_2,
            eccentricity,
            chi_1,
            chi_2,
            luminosity_distance,
            theta_jn,
            phase,
            **kwargs
        )
    else:
        raise RuntimeError("Either TimeDomain or FrequencyDomain should be true")


def native_frequency_domain(
    frequency_array,
    mass_1,
    mass_2,
    eccentricity,
    chi_1,
    chi_2,
    luminosity_distance,
    theta_jn,
    phase,
    **kwargs
):

    if "maximum_frequency" not in kwargs:
        maximum_frequency = frequency_array[-1]
    else:
        maximum_frequency = kwargs["maximum_frequency"]
    if "minimum_frequency" not in kwargs:
        minimum_frequency = frequency_array[0]
    else:
        minimum_frequency = kwargs["minimum_frequency"]

    frequency_bounds = (frequency_array >= minimum_frequency) * (
        frequency_array <= maximum_frequency
    )

    df = frequency_array[1] - frequency_array[0]

    pars = {
        "M": mass_1 + mass_2,
        "q": mass_1 / mass_2,
        "ecc": eccentricity,
        "Lambda1": 0.0,
        "Lambda2": 0.0,
        "chi1": chi_1,
        "chi2": chi_2,
        "domain": 1,  # 1 for TD, 1 for FD
        "arg_out": 0,  # Output hlm/hflm. Default = 0
        # "use_mode_lm": kwargs['use_mode_lm'],  # List of modes to use/output through EOBRunPy
        #  "srate_interp": kwargs['sampling_rate'],  # srate at which to interpolate. Default = 4096.
        "srate_interp": int(max(frequency_array) * 2),
        "use_geometric_units": 0,  # Output quantities in geometric units. Default = 1
        "initial_frequency": kwargs[
            "minimum_frequency"
        ],  # in Hz if use_geometric_units = 0, else in geometric units
        "interp_uniform_grid": 1,  # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        "distance": luminosity_distance,  # Mpc,
        "inclination": theta_jn,
        "output_hpc": 0,
        "df": df,
        **kwargs
    }

    f, hp_real, hp_imag, hc_real, hc_imag = EOBRun_module.EOBRunPy(pars)

    h_plus = np.zeros_like(frequency_array, dtype=np.complex)
    h_cross = np.zeros_like(frequency_array, dtype=np.complex)

    h_plus *= frequency_bounds
    h_cross *= frequency_bounds

    nonzero_mask = (frequency_array >= minimum_frequency) & (
        frequency_array <= maximum_frequency
    )

    h_plus.real[nonzero_mask] = hp_real
    h_plus.imag[nonzero_mask] = hp_imag
    h_cross.real[nonzero_mask] = hc_real
    h_cross.imag[nonzero_mask] = hc_imag

    return dict(plus=h_plus, cross=h_cross)


def fourier_transform_time_domain(
    frequency_array,
    mass_1,
    mass_2,
    eccentricity,
    chi_1,
    chi_2,
    luminosity_distance,
    theta_jn,
    phase,
    **kwargs
):

    # raise RuntimeError(f"Value of {frequency_array}")

    # srate_interp should be determined from

    if "maximum_frequency" not in kwargs:
        maximum_frequency = frequency_array[-1]
    else:
        maximum_frequency = kwargs["maximum_frequency"]
    if "minimum_frequency" not in kwargs:
        minimum_frequency = frequency_array[0]
    else:
        minimum_frequency = kwargs["minimum_frequency"]

    frequency_bounds = (frequency_array >= minimum_frequency) * (
        frequency_array <= maximum_frequency
    )

    df = frequency_array[1] - frequency_array[0]
    srate_interp = int(max(frequency_array) * 2)
    dt = 1.0 / (srate_interp)

    pars = {
        "M": mass_1 + mass_2,
        "q": mass_1 / mass_2,
        "ecc": eccentricity,
        "Lambda1": 0.0,
        "Lambda2": 0.0,
        "chi1": chi_1,
        "chi2": chi_2,
        "domain": 0,  # 0 for TD, 1 for FD
        "arg_out": 0,  # Output hlm/hflm. Default = 0
        # "use_mode_lm": kwargs['use_mode_lm'],  # List of modes to use/output through EOBRunPy
        "srate_interp": srate_interp,
        "use_geometric_units": 0,  # Output quantities in geometric units. Default = 1
        "initial_frequency": kwargs.pop(
            "minimum_frequency"
        ),  # in Hz if use_geometric_units = 0, else in geometric units
        "interp_uniform_grid": 1,  # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        "distance": luminosity_distance,  # Mpc,
        "inclination": theta_jn,
        "output_hpc": 0,
        "dt_interp": dt,
        **kwargs,
    }

    # Since EOB performs an ODE integration, we can't control the actual length of the result, which means the df of the FT won't match with df of frequency_array. Its possible to interpolate but I noticed spurious behaviour at the lower frequency. So we'll fix df by shortening the time domain waveform. (This might not be a good idea...)
    T_max = 1.0 / df
    N_max = int(T_max / dt)

    t, hp, hc = EOBRun_module.EOBRunPy(pars)

    if len(t) > N_max:
        t = t[-N_max:]
        hp = hp[-N_max:]
        hc = hc[-N_max:]

        h_plus = taper_timeseries(
            TimeSeries(hp, delta_t=dt), tapermethod="TAPER_START"
        ).to_frequencyseries()

        h_cross = taper_timeseries(
            TimeSeries(hc, delta_t=dt), tapermethod="TAPER_START"
        ).to_frequencyseries()

    elif len(t) <= N_max:
        deficit = N_max - len(t)
        new_hp = np.zeros(N_max)
        new_hc = np.zeros(N_max)

        new_hp[deficit:N_max] = taper_timeseries(TimeSeries(hp, delta_t=dt), tapermethod="TAPER_START")
        new_hc[deficit:N_max] = taper_timeseries(TimeSeries(hc, delta_t=dt), tapermethod="TAPER_START")
        hp, hc = new_hp, new_hc

        h_plus = TimeSeries(new_hp, delta_t=dt).to_frequencyseries()
        h_cross = TimeSeries(new_hc, delta_t=dt).to_frequencyseries()


    h_plus *= frequency_bounds
    h_cross *= frequency_bounds

    #    nonzero_mask = (frequency_array >= minimum_frequency) \
    #                 & (frequency_array <= maximum_frequency)

    #    h_plus.real[nonzero_mask] = hp_real
    # h_plus.imag[nonzero_mask] = hp_imag
    # h_cross.real[nonzero_mask] = hc_real
    # h_cross.imag[nonzero_mask] = hc_imag

    assert np.all(np.array(h_plus.sample_frequencies) == frequency_array)

    return dict(plus=np.array(h_plus), cross=np.array(h_cross))


def time_domain_eob_resum(
    dt,
    minimum_frequency,
    mass_1,
    mass_2,
    eccentricity,
    luminosity_distance,
    theta_jn,
    phase,
    chi_1=None,
    chi_2=None,
    **kwargs
):
    """
    Generate time domain waveform using TEOBResum with eccentricity.

    All of the named parameters are normal parameters for a BBH. The options to kwargs are passed to the eobresum python interface.

    Possible kwargs:
    "ode_abstol" : Absolute ode tolerance
    "ode_reltol" : Relative ode tolerance
    "use_mode_lm": Which modes to use
    """

    srate_interp = 1.0 / dt
    pars = {
        "M": mass_1 + mass_2,
        "q": mass_1 / mass_2,
        "ecc": eccentricity,
        "Lambda1": 0.0,
        "Lambda2": 0.0,
        "chi1": chi_1 or 0.,
        "chi2": chi_2 or 0.,
        "domain": 0,  # 0 for TD, 1 for FD
        "arg_out": 0,  # Output hlm/hflm. Default = 0
        "srate_interp": srate_interp,
        "use_geometric_units": 0,  # Output quantities in geometric units. Default = 1
        "initial_frequency": minimum_frequency
        ,  # in Hz if use_geometric_units = 0, else in geometric units
        "interp_uniform_grid": 1,  # Interpolate mode by mode on a uniform grid. Default = 0 (no interpolation)
        "distance": luminosity_distance,  # Mpc,
        "inclination": theta_jn,
        "output_hpc": 0,
        "dt_interp": dt,
        **kwargs,
    }

    t, hp, hc = EOBRun_module.EOBRunPy(pars)
    return t, hp, hc
