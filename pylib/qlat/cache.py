class Cache(dict):

    pass

cache = Cache()

def clean_cache(cache = cache):
    for key, val in cache:
        if isinstance(val, Cache):
            clean_cache(val)
        else:
            cache.pop(key)

def mk_cache(key, cache = cache):
    assert not (key in cache)
    cache[key] = Cache()
    return cache[key]

def rm_cache(key, cache = cache):
    assert isinstance(cache[key], Cache)
    cache.pop(key)
