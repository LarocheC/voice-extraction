==========================
Filter bank specifications
==========================

Length of the input signal (compulsory)
---------------------------------------
Before running ``setup``, it is compulsory so fill the ``size`` field with

::

	opts{1}.time.size = length(signal);

It must be a power of 2, which leads to optimally fast Fourier Transforms.
If you have ``K`` signals of the same size ``N``, consider stacking them into a ``KxN`` matrix. This wil automatically vectorize the computation and avoid high-level loop overhead.

Under development is a more general architecture that automates padding to the next power of 2, and adapts to all sizes.


Amount of invariance to translation
-----------------------------------
The integer ``T`` is the amount of invariance to translation that you require. It must also be a power of 2.

A typical value for second-order scattering of audio is ``T=8192``, that is 370 ms at a sample rate of 22 kHz. A smaller ``T`` will not integrate full musical notes or full phonemes ; on the contrary, a bigger ``T`` will blur different notes/phonemes together.

The number of octaves in the filter bank is equal to ``J = log2(T)``.
By default, ``T`` is set equal to ``size`` which means that the corresponding scattering representation ``S`` will be fully translation-invariant.


Quality factor
--------------
The quality factor ``max_Q`` of a band-pass filter is defined as the ratio of its center frequency by its bandwidth. Consequently, for a given center frequency, increasing the quality factor will decrease the bandwidth proportionnally, hence yielding a "sharper" band-pass filter in the frequency domain. This increase in frequency sharpness comes at the cost of increasing the support of the filter in the time domain, which may prevent the representation to distinguish consecutive events.

All the wavelets in a filter bank share the same quality factor: this is why we refer to it as a constant-Q filter bank. Note that this toolbox also allows variable-Q filter banks in order to cope with time support limitations (see section below). This is why the quality factor is ``max_Q``.

Typical values for the first order in audio range from 4 to 16.
Typical values for the second order along time are 1 or 2. 
In the context of multivariable scattering, the value 1 is strongly recommended for any derived variable.

A quality factor of 1, corresponding to the so-called 'dyadic' filter bank, is the default.


Maximum scale
-------------
Note that a potential drawback of the constant-Q filterbank is that the time support of the filters is unbounded at the low frequencies. In audio, it is undesirable that acoustic events more than 100 ms apart fall between the same first-order time bin. To address this issue, this toolbox provides a bound ``max_scale`` that restricts the time support, at the cost of decreasing locally the quality factor.

For instance, for ``max_Q = 12`` and a sample rate of 22 kHz, setting ``max_scale = 2048`` (about 93 ms) will provide constant-Q filters for frequencies above Q/max_scale (about 130 Hz) and constant-bandwidth filters below that limit.
Setting ``max_scale = Inf`` will remove the upper bound on the time support and will guarantee that the quality factor is indeed constant throughout the whole frequency range.

By default, ``max_scale`` is set to ``size``, which means that the time support is only limited by the size of the whole signal.


Number of filters per octave
----------------------------
The integer ``nFilters_per_octave`` specified the rational quantization of the ``gamma`` log-scale variable. In order to cover the whole frequency axis, it is compulsory to have

::

	nFilters_per_octave > max_Q

The number of filters in the filter bank is equal to ``nFilters_per_octave * log2(T)``. Henceforth, note that the computational complexity of the computation is linear in the number of filters per octave of each filter bank.