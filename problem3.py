import random as rand

def checkprime(m):

	if m == 2 or m == 3:
		return True
	if m%2 == 0 or m%3 == 0:
		return False
	x = m-1
	r = 0
	while (x%2 == 0):
		x = x/2
		r = r+1
	d = x

	k = 10

	while (k > 0):
		flag = False
		a = rand.randint(2,m-2)
		v = (a**d)%m
		if ((v == 1)or(v == m-1)):
			flag = True
		else:
			ind = r-1
			while(ind > 0):
				v = (v*v)%m
				if (v == m-1):
					return True
				ind = ind-1

		if flag == False:
			return False
		k =k-1
	return True

def WitnessandLiars(m):
	witness = []
	liars = []
	x = m-1
	r = 0
	while (x%2 == 0):
		x = x/2
		r = r+1
	d = x
	for a in range(2,n-1):
		flag = False
		v = (a**d)%m
		if ((v == 1)or(v == m-1)):
			flag = True
			liars.append(a)
		else:
			ind = r-1
			while(ind > 0):
				v = (v*v)%m
				if (v == m-1):
					flag = True
					liars.append(a)
				ind = ind-1

		if flag == False:
			witness.append(a)
			
	if len(witness) > 0:
		print 'Strong liars:',liars
		print 'Strong witness:',witness
	else:
		print 'Strong liars:',witness
		print 'Strong witness:',liars
def Isprime(number):

	if checkprime(number):
		return True

	return False

n = input()

if Isprime(n):
	print n,'is prime.'
else:
	print n,'is composite.'

m = n+1
while not Isprime(m):
	m = m+1

print m,'is smallest prime greater than',n

WitnessandLiars(n)



