import ctypes
import os

dirname = os.path.dirname(__file__)
if dirname == '':
    dirname = './'
path = os.path.join(dirname, "tf2e.so")
rlib = ctypes.cdll.LoadLibrary(path)

def gen(**params):
    import numpy, pycbc.conversions
    from pycbc.types import FrequencySeries
    from ctypes import c_double, c_void_p
    from pycbc.pnutils import megaparsecs_to_meters

    df = params['delta_f']
    fend = 1000 if 'f_final' not in params else params['f_final']
    if fend == 0:
        fend = 1000
    flow = params['f_lower']
    flen = int(fend / df) + 1

    m1 = params['mass1']
    m2 = params['mass2']
    inc = params['inclination']
    ecc = params['eccentricity']
    if ecc < 1e-5:
        ecc = 1e-5

    lc = params['long_asc_nodes']
    lamc = params['coa_phase']
    dist = params['distance'] / 9.7156118319036E-15

    mchirp = float(pycbc.conversions.mchirp_from_mass1_mass2(m1, m2))
    eta = float(pycbc.conversions.eta_from_mass1_mass2(m1, m2))

    hp = numpy.zeros(flen, dtype=numpy.complex128)
    hc = numpy.zeros(flen, dtype=numpy.complex128)

    f = rlib.generate
    f.argtypes = [c_void_p, c_void_p, c_double, c_double, c_double,
                  c_double, c_double, c_double, c_double, c_double, c_double]

    _ = f(hp.ctypes.data, hc.ctypes.data,
          mchirp, eta, inc, ecc, lamc, lc, dist, fend, df)

    # it appears that the plus / cross data is time inverted
    hp = FrequencySeries(hp.conj(), delta_f=df, epoch=-int(1.0 / df))
    hc = FrequencySeries(hc.conj(), delta_f=df, epoch=-int(1.0 / df))

    kmin = int(flow / df)
    hp[:kmin].clear()
    hc[:kmin].clear()

    return hp, hc

def add_me(**kwds):
    kwds['cpu_fd']['TaylorF2e'] = gen
    kwds['filter_time_lengths']['TaylorF2e'] = kwds['filter_time_lengths']['TaylorF2']
