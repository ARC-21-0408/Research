from specutils.fitting import fit_generic_continuum
cont_norm_spec = spec / fit_generic_continuum(spec)(spec.spectral_axis)
f, ax = plt.subplots()  # doctest: +IGNORE_OUTPUT
ax.step(cont_norm_spec.wavelength, cont_norm_spec.flux)  # doctest: +IGNORE_OUTPUT
ax.set_xlim(654*u.nm, 660*u.nm)  # doctest: +IGNORE_OUTPUT +REMOTE_DATA
