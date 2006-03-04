% Filter design
% Relevant links:
% http://www.stanford.edu/~boyd/magdes.html
% http://www.stanford.edu/class/ee364/lectures/filters.pdf
%
%  fir_chebychev_design.m              - Chebychev design of an FIR filter given a desired H(w)
%  one_over_f_filter.m                 - Design a 1/f spectrum shaping (pink-noise) filter
%  equalizer_design.m                  - Equalizer design example
%  iir_mag_design_bandpass_max_atten.m - Maximize stopband attenuation of a bandpass IIR filter
%  fir_lin_phase_lowpass_max_atten.m   - Maximize stopband attenuation of a linear phase lowpass FIR filter
%  fir_mag_design_lowpass_max_atten.m  - Maximize stopband attenuation of a lowpass FIR filter (magnitude design)
%  iir_mag_design_lowpass_max_atten.m  - Maximize stopband attenuation of a lowpass IIR filter
%  fir_lin_phase_lowpass_min_order.m   - Minimize order of a linear phase lowpass FIR filter
%  fir_mag_design_lowpass_min_order.m  - Minimize order of a lowpass FIR filter (magnitude design)
%  fir_lin_phase_lowpass_min_ripple.m  - Minimize stopband ripple of a linear phase lowpass FIR filter
%  fir_lin_phase_lowpass_min_trans.m   - Minimize transition bandwidth of a linear phase lowpass FIR filter
%  spectral_fact.m                     - Spectral factorization using Kolmogorov 1939 approach.
help Contents
