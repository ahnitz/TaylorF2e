from pycbc.waveform import get_td_waveform
from pycbc.filter import match
hp, hc = get_td_waveform(approximant="TaylorF2e", mass1=10, mass2=10, f_lower=30, delta_t=1.0/4096, eccentricity=0)
hp2, hc = get_td_waveform(approximant="TaylorF2e", mass1=10, mass2=10, f_lower=30, delta_t=1.0/4096, eccentricity=.4)


print match(hp, hp2)
