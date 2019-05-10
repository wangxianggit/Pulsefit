# Pulsefit
A digital pulse decomposition algorithm for studies of short-lived alpha decay in nuclear physics.

A pileup pulse decomposition algorithm here for the decomposition of pileup pulses resulting from successive alpha decays or conversion electron decays of
residuals in a Si detector within short decay time small than 15 microseconds.

mapping: Energy calibration;Time of pulse was determined by fast trapezoidal shaping with short and long shaping parameters for detection of pileup 
pulses with small time difference and small amplitude respectively;Initial amplitude was given by the difference of flat top before and after the leading edge.

spulse:The fit algorithm;Alpha traces of mixed source were used to construct standard pulse response or so called superpulse;Amplitudes of each
alpha traces were normalized by its energy to 1.0 MeV which is obtained from trapezoidal shaping and the standard pulse responses were constructed by
averaging a certain number of alpha traces for each strip;Standard pulse response was fitted to an experimental pulse to derive time and energy;
Two-step fit was applied to improve the energy resolution of pileup pulses with time difference less than 500 ns.
