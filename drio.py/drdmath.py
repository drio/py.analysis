import math

def log_it(l, base=2):
  return map(lambda x: math.log(x, base), l)

def s_to_f(l):
  return map(lambda x: float(x), l)

def average(l):
  to_f = s_to_f(l)
  return sum(to_f) * 1.0 / len(to_f)

def std_dev(l):
  import math
  avg = average(l)
  variance = map(lambda x: (x - avg)**2, s_to_f(l))
  return math.sqrt(average(variance))
