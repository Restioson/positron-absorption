# Code: https://code.activestate.com/recipes/577472-halleys-method-for-solving-equations/
# License: https://github.com/ActiveState/code/tree/master/recipes/Python/577472_Halleys_Method_Solving

# Halley's method for solving f(x)=0
# http://en.wikipedia.org/wiki/Halley%27s_method
# FB - 201011265

import mpmath
import math

h = 0.00000001
eps = 10**(-100000000000000000000)
p = 100

# f(x) to solve
def f(x):
    return mpmath.li(x) - p


# general numerical derivative function
def fp(x, k):
    global h
    if k == 0:
        return f(x)
    else:
        return (fp(x + h, k - 1) - fp(x, k - 1)) / h


# main
def solve_fx(initial):
    x = initial  # initial value

    while True:
        fx = f(x)
        fpx = fp(x, 1)
        '''
        xnew = x - (2.0 * fx * fpx) / (2.0 * fpx * fpx - fx * fp(x, 2))
        '''
        xnew = x - (mpmath.li(x) - p)*(math.log(x))/(1 + (mpmath.li(x) - p)/(2*x))
        if abs(xnew - x) <= eps:
            return xnew

        x = xnew


mpmath.mp.dps = 500

print(math.log(solve_fx(mpmath.mpf(p))))

print(mpmath.ei(math.log(solve_fx(mpmath.mpf(p)))))
'''
print(solve_fx(mpmath.mpf(p)))

print(mpmath.li(solve_fx(mpmath.mpf(p))))
'''
