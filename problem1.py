import hashlib


def hash_short(message, length=16):
    return hashlib.sha1(message).hexdigest()[:length/4]

m1 = "apple sauce for nishna"
h_value = hash_short(m1)
print h_value,m1
f = open('corpus_sherlock.txt','r')
lines = f.readlines()

#print lines
for line in lines:
	#print line
	words = line.split()
	for i in range(len(words)):
		for j in range(i+1,len(words)+1):
			m2 = " ".join(words[i:j])
			#print m2
			h = hash_short(m2)
			if h == h_value:
				print h,m2
