import math

p = 13

h = 8
g = 5

t = int(math.ceil(math.sqrt(p)))

A = []
B = []
for ind in range(t):
	
	A.append(g**(ind*t))
	val = g**ind
	B.append((pow(val, p-2, p)*h)%p)

C = list(set(A)&set(B))

if len(C) > 0:
	i = A.index(C[0])
	j = B.index(C[0])

x = i*t+j

print x





