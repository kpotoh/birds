import numpy as np
# from memory_profiler import profile

# @profile
# def my_func():
#     a = [1] * (10 ** 6)
#     b = [2] * (2 * 10 ** 7)
#     del b
#     return a


def another_func():
    a = np.random.random((1, 120000))
    b = np.random.random((120000, 1347))
    print(a @ b)


def f3():
    return 7 ** 6 * 13 ** 3


if __name__ == '__main__':
    # my_func()
    another_func()
    f3()


# python -m cProfile -s tottime  scripts/check_profile.py | head -n 20