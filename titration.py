import math
from functools import partial

ITER = 10000

def bigFunction(a,b,A,X,w,x):
	return b*(x**4) + (w + b*a + b*X) * (x**3) + (w*A + b*a*X - b*w - a*b*A) * (x**2) - (w**2 + b*a*w + a*A*w) * x - (a*(w**2))

def derivative(a,b,A,X,w,x):
	return 4*b*(x**3) + 3*(w + b*a + b*X) * (x**2) + 2*(w*A + b*a*X - b*w - a*b*A)*x - (w**2 + b*a*w + a*A*w)

def newtonSolver(f, deriv, init):
	x = init
	for i in range(ITER):
		x = x - f(x) / deriv(x)

	return x

def titrationPH(ka, kb, kw, va, vb, ca, cb):
	A = (va * ca) / (va + vb)
	X = (vb * cb) / (va + vb)
	return -math.log10(newtonSolver(partial(bigFunction, ka, kb, A, X, kw), partial(derivative, ka, kb, A, X, kw), 1))

with open("titration.txt", "w") as output:
	for i in range(100):
		output.write(str(0.2*i) + "," + str(titrationPH(1.8e-5, 1e7, 1e-14, i*0.2, 10, 0.1, 0.1)) + "\n")