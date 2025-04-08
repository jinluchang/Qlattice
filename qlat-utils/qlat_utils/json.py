import json
import numpy as np

"""
Mathspp Blog
https://mathspp.com/blog/custom-json-encoder-and-decoder
#
With modifications to add np.ndarray support
"""

class ExtendedEncoder(json.JSONEncoder):
    def default(self, obj):
        name = type(obj).__name__
        try:
            encoder = getattr(self, f"encode_{name}")
        except AttributeError:
            return super().default(obj)
        else:
            encoded = encoder(obj)
            encoded["__extended_json_type__"] = name
            return encoded

class ExtendedDecoder(json.JSONDecoder):
    def __init__(self, **kwargs):
        kwargs["object_hook"] = self.object_hook
        super().__init__(**kwargs)
    def object_hook(self, obj):
        try:
            name = obj["__extended_json_type__"]
            decoder = getattr(self, f"decode_{name}")
        except (KeyError, AttributeError):
            return obj
        else:
            return decoder(obj)

class QlatEncoder(ExtendedEncoder):
    def encode_complex(self, c):
        return {"real": c.real, "imag": c.imag}
    def encode_complex128(self, c):
        return {"real": c.real.item(), "imag": c.imag.item()}
    def encode_complex64(self, c):
        return {"real": c.real.item(), "imag": c.imag.item()}
    def encode_float32(self, f):
        return {"value": f.item()}
    def encode_int64(self, i):
        return {"value": i.item()}
    def encode_int32(self, i):
        return {"value": i.item()}
    def encode_ndarray(self, arr):
        return {"value": arr.tolist()}
    def encode_range(self, r):
        return {"start": r.start, "stop": r.stop, "step": r.step}

class QlatDecoder(ExtendedDecoder):
    def decode_complex(self, obj):
        return complex(obj["real"], obj["imag"])
    def decode_complex128(self, obj):
        return np.complex128(complex(obj["real"], obj["imag"]))
    def decode_complex64(self, obj):
        return np.complex64(complex(obj["real"], obj["imag"]))
    def decode_float32(self, obj):
        return np.float32(obj["value"])
    def decode_int64(self, obj):
        return np.int64(obj["value"])
    def decode_int32(self, obj):
        return np.int32(obj["value"])
    def decode_ndarray(self, obj):
        return np.array(obj["value"])
    def decode_range(self, obj):
        return range(obj["start"], obj["stop"], obj["step"])

def json_dumps(obj, *, indent=None):
    """
    return str dumped from `obj`
    #
    indent (int | str | None)
    """
    return json.dumps(obj, indent=indent, cls=QlatEncoder)

def json_loads(s):
    """
    return obj loaded from str `s`
    """
    return json.loads(s, cls=QlatDecoder)
