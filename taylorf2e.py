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

from pycbc.types import FrequencySeries, zeros, TimeSeries
import numpy

def multi_band(bands=[], lengths=[], startclear=0, **p):
    from pycbc.waveform import get_fd_waveform
    df = p['delta_f']
    fmax = p['f_final']
    flow = p['f_lower']

    bands = [flow] + bands + [fmax]
    dfs = [df] + [1.0 / l for l in lengths]

    dt = 1.0 / (2.0 * fmax)
    tlen = int(1.0 / dt / df)
    flen = tlen / 2 + 1
    wf = TimeSeries(zeros(tlen, dtype=numpy.float32), copy=False, delta_t=dt, epoch=-1.0/df)

    for i in range(len(lengths)+1):
        start = bands[i]
        stop = bands[i+1]
        p2 = p.copy()
        p2['delta_f'] = dfs[i]
        p2['f_lower'] = start
        p2['f_final'] = stop
        hp, hc = get_fd_waveform(**p2)
        hp = hp.astype(numpy.complex64)
        tlen = int(1.0 / dt / dfs[i])
        flen = tlen / 2 + 1
        hp.resize(flen)
        ht = hp.to_timeseries()
        ht[0:int(startclear * ht.sample_rate)] = 0

        wf[len(wf)-len(ht):] += ht

    return wf


def fast_tf2e(**params):

    if 'approximant' in params:
        params.pop('approximant')
    wf = multi_band(bands=[80], lengths=[80], approximant="TaylorF2e", startclear=10, **params)

    maxlen=512
    kmin = int(len(wf) - maxlen * wf.sample_rate)
    if kmin > 0:
        wf[:kmin] = 0

    hp = wf.to_frequencyseries()
    return hp, None

def add_me(**kwds):
    kwds['cpu_fd']['TaylorF2e'] = gen
    kwds['filter_time_lengths']['TaylorF2e'] = kwds['filter_time_lengths']['TaylorF2']

    kwds['cpu_fd']['TaylorF2eB'] = fast_tf2e
    kwds['filter_time_lengths']['TaylorF2eB'] = kwds['filter_time_lengths']['TaylorF2']

