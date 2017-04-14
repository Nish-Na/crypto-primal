import math



p = 13	#Prime number
h = 8	#result of g^xmodp
g = 5	#base

print 'p:',p
print 'h:',h
print 'g:',g

#Calculating the squareroot
t = int(math.ceil(math.sqrt(p)))

A = []
B = []
#finding the values of g^(i*sqrt(p)) and h*g^-j by itirating over all values of i and j
for ind in range(t):
	
	A.append(g**(ind*t))
	val = g**ind
	B.append((pow(val, p-2, p)*h)%p)

#finding the common values
C = list(set(A)&set(B))

if len(C) > 0:
	i = A.index(C[0])
	j = B.index(C[0])

#Calculating value 'x'
x = i*t+j

print x





