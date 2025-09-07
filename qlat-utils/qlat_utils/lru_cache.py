from collections import OrderedDict

class LRUCache:

    """
    Edited from: https://www.geeksforgeeks.org/lru-cache-in-python-using-ordereddict/
    """

    def __init__(self, capacity: int):
        """
        Initializing capacity
        """
        self.cache = OrderedDict()
        self.capacity = capacity

    def clear(self):
        self.cache.clear()

    def __contains__(self, key):
        return key in self.cache

    def get(self, key, default=None):
        if key not in self.cache:
            return default
        else:
            self.cache.move_to_end(key)
            return self.cache[key]

    def __getitem__(self, key):
        """
        We return the value of the key
        that is queried in O(1) and return -1 if we
        don't find the key in out dict / cache.
        And also move the key to the end
        to show that it was recently used.
        """
        if key not in self.cache:
            raise KeyError
        else:
            self.cache.move_to_end(key)
            return self.cache[key]

    def __setitem__(self, key, value) -> None:
        """
        First, we add / update the key by conventional methods.
        And also move the key to the end to show that it was recently used.
        But here we will also check whether the length of our
        ordered dictionary has exceeded our capacity,
        If so we remove the first key (least recently used)
        """
        self.cache[key] = value
        self.cache.move_to_end(key)
        if len(self.cache) > self.capacity:
            self.cache.popitem(last=False)
