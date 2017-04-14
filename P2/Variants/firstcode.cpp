#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

#define PI  3.14159265358979323846264338327950288

struct complx
{
	double real;
	double img;
	complx(double arg1=0,double arg2=0)
	{
		real = arg1,img = arg2;
	}
};

inline complx operator+(const complx& lhs, const complx& rhs)
{
	return complx(lhs.real+rhs.real,lhs.img+rhs.img); 
}

inline complx operator-(const complx& lhs, const complx& rhs)
{
	return complx(lhs.real-rhs.real,lhs.img-rhs.img); 
}

inline complx operator*(const complx& lhs, const complx& rhs)
{
	return complx(lhs.real*rhs.real - lhs.img*rhs.img, lhs.real*rhs.img + lhs.img*rhs.real); 
}

inline complx operator*(const double&scale, const complx& val)
{
	return complx(scale*val.real, scale*val.img); 
}

inline complx operator*(const complx& val,const double& scale)
{
	return complx(scale*val.real, scale*val.img); 
}

inline complx operator/(const complx& val,const double& scale)
{
	return complx(val.real/scale, val.img/scale); 
}

inline complx operator/(const complx& lhs, const complx& rhs)
{
	double scale = (rhs.real*rhs.real)+(rhs.img*rhs.img);
	return complx(lhs.real*rhs.real + lhs.img*rhs.img, lhs.img*rhs.real - lhs.real*rhs.img)/scale; 
}

inline complx omegaval(ll N,ll i)
{
	double rad = (2*PI*i)/(N);
	return complx(cos(rad),sin(rad));
}

#define FOR(i, a, b) for (ll i = (a); i < (b); ++i)
#define REP(i, n) FOR(i, 0, n)
#define LIM 1000000
#define MAX 1100000

/***************** Polynomial functions ******************/
struct poly
{
	ll degree;
	ll* coeff;
	poly(ll deg)
	{
		degree = deg;
		coeff = new ll[deg+1];
	}
};

poly* lhs;
poly* rhs;
poly* mulres;
poly* powres;

void polyprint(poly* P)
{
	// prints polynomial
	ll deg = P->degree;
	for(ll i=deg;i>0;i--)
		cout << P->coeff[i] << "x^" << i << " + ";
	cout << P->coeff[0] << endl;
}

void polymod(poly* P,ll r,ll n)
{
	// stores P mod(x^r-1,n) in P
	ll deg = P->degree;
	while(deg>=r)
	{
		P->coeff[deg-r] += P->coeff[deg];
		if(P->coeff[deg-r]>=n)
			P->coeff[deg-r]-=n;
		P->coeff[deg] = 0;
		deg--;
	}
	while(deg>0 && P->coeff[deg]==0)
		deg--;
	P->degree = deg;
}

ll N;
complx omega[LIM];
complx a1[MAX], a2[MAX];
complx z1[MAX], z2[MAX];

ll A0[MAX], A1[MAX];
ll B0[MAX], B1[MAX];
ll C0[MAX], C1[MAX], C2[MAX];


void fft(complx *a, complx *z, ll m = N)
{
	if(m == 1)
		z[0] = a[0];
	else
	{
		ll s = N / m;
		m /= 2;

		fft(a, z, m);
		fft(a + s, z + m, m);

		REP(i,m)
		{
			complx c = omega[s*i] * z[m + i];
			z[m + i] = z[i] - c;
			z[i] = z[i] + c;
		}
	}
}

void mult(ll *a, ll *b, ll *c, int len)
{
	N = 2 * len;
	while(N & (N - 1))
		++N;
	REP(i, N) a1[i] = complx(0, 0);
	REP(i, N) a2[i] = complx(0, 0);
	REP(i, len) a1[i] = complx(a[i], 0);
	REP(i, len) a2[i] = complx(b[i], 0);

	REP(i, N) omega[i] = omegaval(N,i);
	fft(a1, z1, N);
	fft(a2, z2, N);

	REP(i, N) omega[i] = complx(1, 0)/omega[i];
	REP(i, N) a1[i] = z1[i] * z2[i] / complx(N,0);
	fft(a1, z1, N);

	REP(i, 2 * len) c[i] = round(z1[i].real);
}

void polymulmod(poly* a, poly* b, ll r, ll n)
{
	// n is modulo
	// len stores maximum degree of a,b
	ll len = a->degree;
	if(b->degree>len)
		len = b->degree;
	
	// pad extra zeroes to smaller polynomial
	for(ll i=a->degree+1;i<=len;i++)
		a->coeff[i]=0;
	for(ll i=b->degree+1;i<=len;i++)
		b->coeff[i]=0;

	// make len store number of terms
	len++;

	REP(i, len) A0[i] = a->coeff[i] & 0xFFFF;
	REP(i, len) A1[i] = a->coeff[i] >> 16;

	REP(i, len) B0[i] = b->coeff[i] & 0xFFFF;
	REP(i, len) B1[i] = b->coeff[i] >> 16;

	mult(A0, B0, C0, len);
	mult(A1, B1, C2, len);

	REP(i, len) A0[i] += A1[i];
	REP(i, len) B0[i] += B1[i];
	mult(A0, B0, C1, len);
	REP(i, 2 * len) C1[i] -= C0[i] + C2[i];

	REP(i, 2 * len) C1[i] %= n;
	REP(i, 2 * len) C2[i] %= n;
	REP(i, 2 * len) a->coeff[i] = (C0[i] + (C1[i] << 16) + (C2[i] << 32)) % n;

	// remove extra zeroes and find degree of resulting polynomial
	a->degree = 2*len;
	while(a->degree>0 && a->coeff[a->degree]==0)
		a->degree--;
	
	polymod(a,r,n);
}

/*
void polymulmod(poly* P,poly* Q,ll r,ll n)
{
	// P*Q modulo x^r-1,n
	// stores result in P 
	ll degP = P->degree;
	ll degQ = Q->degree;
	mulres->degree = degP+degQ;
	for(ll i=0;i<=(2*r);i++)
		mulres->coeff[i] = 0;
	for(ll i=0;i<=degP;i++)
	{
		for(ll j=0;j<=degQ;j++)
		{
			mulres->coeff[i+j] += (P->coeff[i]*Q->coeff[j])%n;
			if(mulres->coeff[i+j]>=n)
				mulres->coeff[i+j]-=n;
		}
	}
	polymod(mulres,r,n);
	P->degree = mulres->degree;
	for(ll i=0;i<=mulres->degree;i++)
		P->coeff[i] = mulres->coeff[i];
}
*/

void polypow(poly *P,ll r,ll n)
{
	// stores p**n mod(x^r-1,n) in p
	ll b=n;
	// Initialize powres as 1 to store result of exponentiation
	powres->degree = 0;
	powres->coeff[0] = 1;
	while(b!=0)
	{
		if(b%2==1)
			polymulmod(powres,P,r,n);
		polymulmod(P,P,r,n);
		b/=2;
	}
	P->degree = powres->degree;
	for(ll i=0;i<=powres->degree;i++)
		P->coeff[i] = powres->coeff[i];
	for(ll i=powres->degree+1;i<=(2*r);i++)
		P->coeff[i] = 0;
}

/******************** Integer functions **************************/

ll gcd(ll a, ll b)
{
	ll t;
	while(b!=0)
	{
		t = a;
		a = b;
		b = t%b;
	}
	return a;
}

ll expo(ll a,ll b)
{
	ll res = 1;
	while(b!=0)
	{
		if(b%2==1)
			res*=a;
		a*=a;
		b/=2;
	}
	return res;
}

ll powmod(ll a,ll b,ll m)
{
	ll res = 1;
	while(b!=0)
	{
		if(b%2==1)
			res=(res*a)%m;
		a=(a*a)%m;
		b/=2;
	}
	return res;
}

ll logf(ll N)
{
	assert (N>0);
	ll logn=0;
	while(N>0)
	{
		N/=2;
		logn++;
	}
	return logn-1;
}

ll logc(ll N)
{
	assert (N>0);
	ll logn=0;
	ll bits=0;
	while(N>0)
	{
		if(N%2==1)
			bits++;
		N/=2;
		logn++;
	}
	return logn-(bits==1);
}

ll phi(ll N)
{
	ll result = N;
	for(ll p=2;p*p<=N;p++)
	{
		if(N%p==0)
		{
			while(N%p==0)
				N/=p;
			result-=result/p;
		}
	}
	if(N>1)
		result-=result/N;
	return result;
}

bool IsPerfectPower(ll N)
{
	ll logn,t,left,right,mid,exp;
	logn = logf(N);
	for(t=2;t<=logn;t++)
	{
		left = 2, right = N;
		while(left<=right)
		{
			mid = (left+right)/2;
			exp = expo(mid,t);
			if(exp<N)
				left = mid+1;
			else if(exp>N)
				right = mid-1;
			else
				return true;
		}
	}
	return false;
}

/********************** Main Algorithm ******************************/
bool isequal(poly* P,poly* Q)
{
	// returns true in P and Q are equal
	ll degP = P->degree;
	ll degQ = Q->degree;
	for(;degQ>degP;degQ--)
		if(Q->coeff[degQ]!=0)
			return false;
	for(;degP>degQ;degP--)
		if(P->coeff[degP]!=0)
			return false;
	for(;degP>=0 && P->coeff[degP]==Q->coeff[degP];degP--);
	return (degP==-1);
}
bool checkAKS(ll a,ll n,ll r)
{
	// returns true if (x+a)^n = x^n+a mod(x^r-1,n)
	// polynomial is stored as a + bx + cx^2 + ......

	// Initialize lhs as x+a and compute (x+a)^n mod(x^r-1,n)
	lhs->degree = 1;
	lhs->coeff[1] = 1;
	lhs->coeff[0] = a;
	polypow(lhs,r,n);

	// Initialize rhs as x and compute x^n mod(x^r-1,n) add a and compute x^n+a mod(x^r-1,n)
	rhs->degree = 1;
	rhs->coeff[1] = 1;
	rhs->coeff[0] = 0;
	polypow(rhs,r,n);
	rhs->coeff[0] += (a%n);
	if(rhs->coeff[0] >=n )
		rhs->coeff[0] -= n;
	polymod(rhs,r,n);
	/*
	   cout << "FOR " << a << endl;
	   polyprint(lhs);
	   polyprint(&rhs);
	*/
	return isequal(lhs,rhs);
}

bool AKS(ll N)
{
	// Check if N is a perfect power
	if(IsPerfectPower(N))
		return false;

	// Find smallest r such that ordr(N) >= (logn)^2
	double t = log2((double)N);
	t*=t;
	ll maxk = ceil(t);
	cout << "maxk = " << maxk << endl;
	bool nextR = true;
	ll r,k,temp,a,gcdan;
	for(r=2;nextR;r++)
	{
		if(gcd(r,N)!=1)
			continue;
		nextR = false;
		temp=1;
		for(k=1;(!nextR) && (k<=maxk);k++)
		{
			temp = (temp*N)%r;
			nextR = (temp==1 || temp==0);
		}
	}
	r--;
	cout << "r = " << r << endl;

	// If 1 < gcd(a,n) < n for some a<=r output composite
	ll maxr = r;
	if(r>=N)
		maxr = N-1;
	for(a=maxr;a>1;a--)
	{
		gcdan = gcd(a,N);
		if(gcdan!=1)
			return false;
	}

	// If N<=r output prime
	if(N<=r)
		return true;

	// Find Max a to perform main step of AKS
	ll maxa = ceil((double)sqrt(phi(r))*(double)log2(N));
	cout << "maxa = " << maxa << endl;

	/*******************************************/
	lhs = new poly(2*r);
	lhs->degree = 0;
	rhs = new poly(2*r);
	rhs->degree = 0;
	powres = new poly(2*r);
	powres->degree = 0;
	mulres = new poly(2*r);
	mulres->degree = 0;
	for(ll a=1;a<=maxa;a++)
	{
		cout << a << endl;
		if(!checkAKS(a,N,r))
			return false;
	}
	return true;
}

int main()
{
	ll n;
	n = 99929;
	n = 7069;
	n = 1000000007;
	if(AKS(n))
		cout << n << " is prime." << endl;
	else
		cout << n << " is not prime." << endl;
	return 0;
}
