import cProfile
import pstats
import sys

import numpy as np
# from memory_profiler import profile

# @profile
# def my_func():
#     a = [1] * (10 ** 6)
#     b = [2] * (2 * 10 ** 7)
#     del b
#     return a


def profiler(nlines=10):
    def decorator_profiler(func):
        def wrapper_profiler(*args, **kwargs):
            profile = cProfile.Profile()
            profile.enable()
            value = func(*args, **kwargs)
            profile.disable()
            ps = pstats.Stats(profile, stream=sys.stderr)
            print(f"{'-' * 30}\nFunction: {func.__name__}\n{'-' * 30}", file=sys.stderr)
            ps.sort_stats('cumtime', 'calls')
            ps.print_stats(nlines)
            return value
        return wrapper_profiler
    return decorator_profiler


@profiler()
def another_func():
    a = np.random.random((1, 120000))
    b = np.random.random((120000, 1347))
    print(a @ b)


@profiler(3)
def f3():
    return 7 ** 6 * 13 ** 3


if __name__ == '__main__':
    # my_func()
    another_func()
    f3()


# python -m cProfile -s tottime  scripts/check_profile.py | head -n 20