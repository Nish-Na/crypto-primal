import random

p = 37 	# prime number
g = 5	# Base

a = random.randint(2,1000000) 	# secret key of Alice 
A = (g**a)%p;	# generating public key of Alice

b = random.randint(2,1000000) # secret key of Bob 
B = (g**b)%p; 	# generating public key of Bob


print 'public key of Alice is ',A,', public key of Bob is ',B

# session keys
s1 = (B**a)%p 
s2 = (A**b)%p

print 'The session keys are s1=',s1,', s2=',s2

