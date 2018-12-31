import numpy
from matplotlib import pyplot
import sympy

sympy.init_printing()

x = sympy.symbols('x')

func = (sympy.cos(x)**2 * sympy.sin(x)**3)/(4*x**5 * sympy.exp(x))

funcprime = func.diff(x)

funcderv = sympy.lambdify((x), funcprime)

print(funcderv(2.2))