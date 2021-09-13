import cqlat as c

from qlat.field import *

import math

def mk_phase_field(geo: Geometry, lmom):
    # lmom is in lattice momentum unit
    # exp(i * 2*pi/L * lmom \cdot xg )
    f = Field("Complex", geo, 1)
    c.set_phase_field(f, lmom)
    return f

class FastFourierTransform:

    def __init__(self, fft_infos, *, is_normalizing = False):
        # fft_infos = [ ( fft_dir, is_forward, ), ... ]
        self.fft_infos = fft_infos
        self.is_normalizing = is_normalizing

    def __mul__(self, field):
        for fft_dir, is_forward in self.fft_infos:
            f = field.copy(False)
            c.fft_dir_complex_field(f, field, fft_dir, is_forward)
            field = f
        if self.is_normalizing and self.fft_infos:
            total_site = field.total_site()
            scale_factor = 1
            for fft_dir, is_forward in self.fft_infos:
                scale_factor *= total_site[fft_dir]
            scale_factor = 1.0 / math.sqrt(scale_factor)
            field *= scale_factor
        return field

def mk_fft(is_forward, *, is_only_spatial = False, is_normalizing = False):
    if is_only_spatial:
        fft_infos = [
                (0, is_forward,),
                (1, is_forward,),
                (2, is_forward,),
                ]
        return FastFourierTransform(fft_infos, is_normalizing = is_normalizing)
    else:
        fft_infos = [
                (0, is_forward,),
                (1, is_forward,),
                (2, is_forward,),
                (3, is_forward,),
                ]
        return FastFourierTransform(fft_infos, is_normalizing = is_normalizing)
