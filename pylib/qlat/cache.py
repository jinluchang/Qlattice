class Cache(dict):

    pass

cache = Cache()

def clean_cache(cache = cache):
    for key, val in list(cache.items()):
        if isinstance(val, Cache):
            clean_cache(val)
        else:
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
