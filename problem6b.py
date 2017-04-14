import random

p = 37 	# prime number
g = 5	# Base


a = random.randint(2,1000000) 	# secret key of Alice 
A = (g**a)%p;	# generating public key of Alice

b = random.randint(2,1000000) # secret key of Bob 
B = (g**b)%p; 	# generating public key of Bob

c = random.randint(2,1000000) # secret key of MIM 
C = (g**c)%p; 	# generating public key of MIM


print 'public key of Alice is ',A,', public key of MIM is ',C,', public key of Bob is ',B

# session keys A-C
s1 = (C**a)%p 
s2 = (A**c)%p

print 'The session keys of A&C are s1=',s1,', s2=',s2

# session keys B-C
s3 = (C**b)%p 
s4 = (B**c)%p

print 'The session keys of B&C are s3=',s3,', s4=',s4