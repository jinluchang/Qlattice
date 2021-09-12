from qlat.timer import *

class Cache(dict):

    def __init__(self, *keys):
        self.keys = keys

cache = Cache()

@timer
def list_cache(cache = cache):
    l = dict()
    for key, val in cache.items():
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
def clean_cache(cache = cache):
    info_str = show_cache_keys(cache.keys)
    displayln_info(f"clean_cache: cache{info_str}:")
    for key, val in list(cache.items()):
        if isinstance(val, Cache):
            clean_cache(val)
        else:
            displayln_info(f"clean_cache: cache{info_str}['{key}']")
            cache.pop(key)

def mk_cache(*keys, cache = cache):
    # if not exist, otherwise return existing elements
    assert keys
    for key in keys:
        if key in cache:
            cache = cache[key]
            assert isinstance(cache, Cache)
        else:
            cache[key] = Cache(*(cache.keys), key)
            cache = cache[key]
    return cache

@timer
def rm_cache(*keys, cache = cache):
    # if exist
    assert keys
    for key in keys[-1]:
        if key not in cache:
            return
        cache = cache[key]
        assert isinstance(cache, Cache)
    key = keys[-1]
    if key not in cache:
        return
    assert isinstance(cache[key], Cache)
    cache.pop(key)
