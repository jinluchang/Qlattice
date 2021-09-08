from qlat.timer import *

class Cache(dict):

    pass

cache = Cache()

@timer
def list_cache(cache = cache):
    l = []
    for key, val in cache.items():
        if isinstance(val, Cache):
            l.append(( key, list_cache(val), ))
        else:
            l.append(key)
    return l

@timer
def clean_cache(cache = cache, info_tag = ""):
    displayln_info(f"clean_cache: {info_tag}:")
    for key, val in list(cache.items()):
        if isinstance(val, Cache):
            clean_cache(val, f"{info_tag}/'{key}'")
        else:
            displayln_info(f"clean_cache: {info_tag}/'{key}'")
            cache.pop(key)

def mk_cache(*keys, cache = cache):
    # if not exist, otherwise return existing elements
    assert keys
    for key in keys[:-1]:
        if key in cache:
            cache = cache[key]
        else:
            cache[key] = Cache()
            cache = cache[key]
    key = keys[-1]
    if key not in cache:
        cache[key] = Cache()
    else:
        assert isinstance(cache[key], Cache)
    return cache[key]

@timer
def rm_cache(*keys, cache = cache):
    # if exist
    assert keys
    for key in keys[-1]:
        if key not in cache:
            return
        cache = cache[key]
    key = keys[-1]
    if key not in cache:
        return
    assert isinstance(cache[key], Cache)
    cache.pop(key)
