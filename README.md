# eobresum_bilby

Thin wrapper around the 'eccentric' branch of https://bitbucket.org/eob_ihes/teobresums/src/master/ which allows compatibility with Bilby https://git.ligo.org/lscsoft/bilby. 

The functions defined are: 

`eccentric_binary_black_hole_eob_resum_nonspinning` satisfies the interface for a `frequency_domain_source_model` used by Bilby's waveform generators. Takes extra keyword arguments of `'TimeDomain'` or `'FrequencyDomain'`, and any keyword arguments used by EOBResum, e.g. `ode_abstol, ode_reltol, use_mode_lm`

`eccentric_binary_black_hole_eob_resum_aligned_spin` The same function which allows aligned spins. 

`native_frequency_domain` Uses the frequency domain model provided by EOBResum NOT RECOMMENDED, doesn't seem to produce correct waveforms. 

`fourier_transform_time_domain` Uses time domain model provided by EOBResum and fourier transforms the result before returning to Bilby. 

`time_domain_eob_resum`: Produces time domain EOBResum waveforms. 

Example: 

``` python 
kwargs = {'ode_abstol': 1e-10, 
          'ode_reltol': 1e-9, 
          'use_mode_lm': [1] #(2,2) mode
          }

hp, hx = time_domain_eob_resum(dt = 2048., minimum_frequency=15, mass_1=30., mass_2=30, eccentricity=0.1, luminosity_distance=500, theta_jn = 0., phase=0., 
                              chi_1 = 0.1, chi_2 = -0.2, **kwargs)
