from qlat_utils.timer import *

import qlat_utils.c as c

from qlat_utils.c import get_all_caches_info

class Cache(dict):

    def __init__(self, *keys):
        super().__init__()
        self.keys = keys

###

cache = Cache()

@timer
def list_cache(ca = cache):
    l = dict()
    for key, val in ca.items():
        if isinstance(val, Cache):
            l[key] = list_cache(val)
        else:
            l[key] = None
    return l

def show_cache_keys(keys):
    if not keys:
        return ""
    else:
        return "['" + "']['".join(keys) + "']"

@timer
def clean_cache(ca = cache):
    """
    Remove values of cache, but keep all the structures
    """
    info_str = show_cache_keys(ca.keys)
    items = list(ca.items())
    displayln_info(0, f"clean_cache: cache{info_str}: len={len(items)}")
    for key, val in items:
        if isinstance(val, Cache):
            clean_cache(val)
        else:
            # displayln_info(1, f"clean_cache: cache{info_str}['{key}']")
            ca.pop(key)

def mk_cache(*keys, ca = cache):
    """
    make cache if it does not exist, otherwise return existing elements
    """
    assert keys
    for key in keys:
        if key in ca:
            ca = ca[key]
            assert isinstance(ca, Cache)
        else:
            ca[key] = Cache(*(ca.keys), key)
            ca = ca[key]
    return ca

@timer
def rm_cache(*keys, ca = cache):
    """
    remove cache if it exist
    """
    assert keys
    for key in keys[-1]:
        if key not in ca:
            return
        ca = ca[key]
        assert isinstance(ca, Cache)
    key = keys[-1]
    if key not in ca:
        return
    assert isinstance(ca[key], Cache)
    ca.pop(key)

def get_all_caches_info():
    return c.get_all_caches_info()

@timer
def clear_all_caches():
    """
    clean python level cache and then C++ level cache
    """
    clean_cache()
    c.clear_all_caches()
