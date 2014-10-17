import math

def log_it(l, base=2):
    return map(lambda x: math.log(x, base), l)

def s_to_f(l):
    return map(lambda x: float(x), l)

def average(l):
    to_f = s_to_f(l)
    return sum(to_f) * 1.0 / len(to_f)

def std_dev(l):
    avg = average(l)
    variance = map(lambda x: (x - avg)**2, s_to_f(l))
    return math.sqrt(average(variance))

def abs_val(n):
    if n >= 0:
        return n
    return -1 * n

def is_prime(n):
    if n == 0:
        return False
    if n == 1:
        return True
    if n % 2 == 0 and n > 2:
        return False
    for i in range(3, int(math.sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True
