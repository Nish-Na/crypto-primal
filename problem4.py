import random
import hashlib
import string


def key_generator():
	return ''.join(random.choice(string.ascii_uppercase + string.digits + string.ascii_lowercase) for _ in range(256))

def commit(msg,key):
	return hashlib.sha256(msg+key).hexdigest()

def verify(com,key,msg):
	if(hashlib.sha256(msg+key).hexdigest() == com):
		print "verified"
	else:
		print "Not verified"

def gcd(p, q):
    while q != 0:
        p, q = q, p % q
    return p

def modmulinv(e,z):
	d = 0
	x1 = 0
	x2 = 1
	y1 = 1
	temp_b = z

	while e > 0:
	    t1 = temp_b/e
	    t2 = temp_b - t1 * e
	    temp_b = e
	    e = t2
	    
	    x = x2- t1* x1
	    y = d - t1 * y1
	    
	    x2 = x1
	    x1 = x
	    d = y1
	    y1 = y

	if temp_b == 1:
	    return d + z 

def rsa_generatekey():
	p=9769
	q=9781
	n = p*q
	z = (p-1)*(q-1)
	g = 0
	while g!=1:
		e = random.randrange(1,z)
		g = gcd(z,e)
	d = modmulinv(e, z)
	return ((e,n),(d,n))

def rsa_encrypt(com,key):
	e, n = key
	cipher = [pow(ord(x),e,n) for x in com]
	return cipher

def rsa_verify(cipher,com,key):
	d, n = key
	# print d,n
	# print cipher
	plain = [chr(pow(num,d,n)) for num in cipher]
	# print com
	if (com == ''.join(plain)):
		print "RSA-Verified"
	else:
		print "RSA-not Verified"



def main():
	#message
	msg = "Dont go away"
	#random key
	print 'message is:',msg
	print '\nGenerating key ....'
	key = key_generator()
	print "key is:",key
	#print "msg is",msg
	#committing message
	print "\ncommiting message ...."
	com = commit(msg,key)
	print "Commited"
	print 'com:',com
	#verifying message
	print "\nVerifying Commit ...."
	verify(com,key,msg)
	print '\n----------------------------------------------------------'
	#RSA-encryption
	#key_pair
	print "\nRSA_Encryption"
	print "\nGenerating key pair ...."
	sk,pk = rsa_generatekey()
	print 'private key is',sk,' Public key is',pk
	#encryption
	print "\nEncrypting ...."
	cipher = rsa_encrypt(com,sk)
	print "Encrypted"
	print 'cipher:',cipher
	#Decryption & Verification
	print "\nSignature Verifying ...."
	rsa_verify(cipher,com,pk)


if __name__ == '__main__':
	main()