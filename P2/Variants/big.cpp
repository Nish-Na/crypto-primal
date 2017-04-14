#include<bits/stdc++.h>
#include<ttmath/ttmath.h>
using namespace std;
typedef long long int ll;
typedef unsigned long long int ull;

typedef ttmath::Int<10> bigint;
typedef ttmath::Big<10,5> bigdouble;

#define PI  3.14159265358979323846264338327950288
struct complx
{
	bigdouble real;
	bigdouble img;
	complx(bigdouble arg1=0,bigdouble arg2=0)
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

inline complx operator*(const bigdouble&scale, const complx& val)
{
	return complx(scale*val.real, scale*val.img); 
}

inline complx operator*(const complx& val,const bigdouble& scale)
{
	return complx(scale*val.real, scale*val.img); 
}

inline complx operator/(const complx& val,const bigdouble& scale)
{
	return complx(val.real/scale, val.img/scale); 
}

inline complx operator/(const complx& lhs, const complx& rhs)
{
	bigdouble scale = (rhs.real*rhs.real)+(rhs.img*rhs.img);
	return complx(lhs.real*rhs.real + lhs.img*rhs.img, lhs.img*rhs.real - lhs.real*rhs.img)/scale; 
}

inline complx omegaval(ll N,ll i)
{
	bigdouble rad = (double)(2.0*PI*i)/(N);
	return complx(ttmath::Cos(rad),ttmath::Sin(rad));
}


void lltobig(ll& lhs,bigint& rhs)
{
	// assign long long to bigint by converting to str in middle
	stringstream ss;
	ss << lhs;
	rhs = ss.str();
}

void bigtoll(bigint& lhs,ll& rhs)
{
	// assign bigint to long long int by converting to str in middle
	stringstream ss(lhs.ToString());
	ss >> rhs;
}

double lognd; // stores logn precisely
int lognf; // stores floor of logn
int lognc; // stores ceil of logn

void init(bigint& N)
{
	bigdouble temp;
	temp.Log(N,2);
	lognd = temp.ToDouble();
	lognf = floor(lognd);
	lognc = ceil(lognd);
	cout << "precise logn = " << lognd << endl;
	cout << "floor of logn = " << lognf << endl;
	cout << "ceil of logn = " << lognc << endl;
}

ll gcd(bigint N, bigint r)
{
	// using the fact gcd(a,b) = gcd(b,a%b) create a,b from N,r so that no need to use bigint
	ll t,a,b;
	N = N%(r.ToInt());
	a = r.ToInt(),b = N.ToInt();
	while(b!=0)
	{
		t = a;
		a = b;
		b = t%b;
	}
	return a;
}

bigint expo(bigint a,int b)
{
	bigint res = 1;
	while(b!=0)
	{
		if(b&1)
			res*=a;
		a*=a;
		b>>=1;
	}
	return res;
}

ll phi(ll N)
{
	ll result = N;
	ll p;
	for(p=2;p*p<=N;p++)
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

bool IsPerfectPower(bigint& N)
{
	int t;
	bigint left,right,mid,exp;
	for(t=2;t<=lognf;t++)
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

/******************************************************************************************************************/
#define FOR(i, a, b) for (ll i = (a); i < (b); ++i)
#define REP(i, n) FOR(i, 0, n)
#define LIM 1000000
#define MAX 1100000
struct poly
{
	ll degree;
	bigint* coeff;
	poly(ll deg)
	{
		degree = deg;
		coeff = new bigint[deg+1];
	}
};

poly* lhs;
poly* rhs;
poly* mulres;
poly* powres;

void polyprint(poly* P)
{
	ll deg = P->degree;
	for(ll i=deg;i>0;i--)
		cout << P->coeff[i] << "x^" << i << " + ";
	cout << P->coeff[0] << endl;
}

void polymod(poly* P,ll r,bigint n)
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

ll Nfft;
complx omega[LIM];
complx a1[MAX], a2[MAX];
complx z1[MAX], z2[MAX];

bigint A0[MAX], A1[MAX];
bigint B0[MAX], B1[MAX];
bigint C0[MAX], C1[MAX], C2[MAX];


void fft(complx *a, complx *z, ll m = Nfft)
{
	if(m == 1)
		z[0] = a[0];
	else
	{
		ll s = Nfft / m;
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

void mult(bigint *a, bigint *b, bigint *c, ll len)
{
	Nfft = 2 * len;
	while(Nfft & (Nfft - 1))
		++Nfft;

	REP(i, Nfft) a1[i] = complx(0, 0);
	REP(i, Nfft) a2[i] = complx(0, 0);
	REP(i, len) a1[i] = complx(a[i], 0);
	REP(i, len) a2[i] = complx(b[i], 0);

	REP(i, Nfft) omega[i] = omegaval(Nfft,i);
	fft(a1, z1, Nfft);
	fft(a2, z2, Nfft);

	REP(i, Nfft) omega[i] = complx(1, 0)/omega[i];
	REP(i, Nfft) a1[i] = z1[i] * z2[i] / complx((double)Nfft,0);
	fft(a1, z1, Nfft);

	bigdouble temp;
	REP(i, 2 * len)
	{
		temp = ttmath::Round((z1[i].real));
		c[i] = temp.ToString();
	}
}

/*
void polymulmod(poly* a, poly* b, ll r, bigint& n)
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

	REP(i, len)
	{
		A0[i] = a->coeff[i];
		A0[i].BitAnd(0xFFFF);
	}
	REP(i, len)
	{
		A1[i] = a->coeff[i];
		A1[i].Rcr(16,0);
	}

	REP(i, len)
	{
		B0[i] = b->coeff[i];
		B0[i].BitAnd(0xFFFF);
	}
	REP(i, len)
	{
		B1[i] = b->coeff[i];
		B1[i].Rcr(16,0);
	}

	mult(A0, B0, C0, len);
	mult(A1, B1, C2, len);

	REP(i, len) A0[i] += A1[i];
	REP(i, len) B0[i] += B1[i];
	mult(A0, B0, C1, len);
	REP(i, 2 * len) C1[i] -= C0[i] + C2[i];

	REP(i, 2 * len) C1[i] %= n;
	REP(i, 2 * len) C2[i] %= n;
	REP(i, 2 * len)
	{
		C1[i].Rcl(16,0);
		C2[i].Rcl(32,0);
		a->coeff[i] = (C0[i] + C1[i] + C2[i]) % n;
	}

	// remove extra zeroes and find degree of resulting polynomial
	a->degree = 2*len;
	while(a->degree>0 && a->coeff[a->degree]==0)
		a->degree--;

	polymod(a,r,n);
}
*/

void polymulmod(poly* P,poly* Q,ll r,bigint& n)
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


void polypow(poly *P,ll r,bigint& n)
{
	// stores p**n mod(x^r-1,n) in p
	bigint b = n;
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

bool checkAKS(ll a,bigint& n,ll r)
{
	// returns true if (x+a)^n = x^n+a mod(x^r-1,n)
	// polynomial is stored as a + bx + cx^2 + ......

	// Initialize lhs as x+a and compute (x+a)^n mod(x^r-1,n)
	bigint biga;
	lltobig(a,biga);
	lhs->degree = 1ll;
	lhs->coeff[1] = 1;
	lhs->coeff[0] = biga;
	polypow(lhs,r,n);

	// Initialize rhs as x and compute x^n mod(x^r-1,n) add a and compute x^n+a mod(x^r-1,n)
	rhs->degree = 1;
	rhs->coeff[1] = 1;
	rhs->coeff[0] = 0;
	polypow(rhs,r,n);
	rhs->coeff[0] += (biga%n);
	if(rhs->coeff[0] >=n )
		rhs->coeff[0] -= n;
	polymod(rhs,r,n);
	/*
	cout << "FOR " << a << endl;
	polyprint(lhs);
	polyprint(rhs);
	*/
	return isequal(lhs,rhs);
}

/******************************************************************************************************************/
bool AKS(bigint& N)
{
	// Check if N is a perfect power
	if(IsPerfectPower(N))
		return false;

	// Find smallest r such that ordr(N) >= (logn)^2
	ll maxk = ceil(lognd*lognd);
	bool nextR = true;
	bigint N2,bigr;
	ll Nmodr,r,temp,k;

	for(bigr=2,r=2;nextR;bigr++,r++)
	{
		if(gcd(N,bigr)!=1)
			continue;
		nextR = false;
		temp = 1;
		N2 = N%bigr;
		bigtoll(N2,Nmodr); // stores N mod r to avoid bigint at next steps
		for(k=1;(!nextR) && (k<=maxk);k++)
		{
			temp = (temp*Nmodr)%r;
			nextR = (temp==1 || temp==0);
		}
	}
	r--;

	// If 1 < gcd(a,n) < n for some a<=r output composite
	ll maxr, gcdan;
	bigtoll(bigr,maxr);
	if(bigr>N)
	{
		bigtoll(N,maxr);
		maxr--; // stores maxr = N-1
	}
	lltobig(maxr,N2);
	for(;N2>1;N2--)
	{
		gcdan = gcd(N,N2);
		if(gcdan!=1)
			return false;
	}

	// If N<=r output prime
	if(N<=bigr)
		return true;

	// Find Max a to perform main step of AKS
	ll maxa = ceil(sqrtl(phi(r))*lognd);
	cout << "maxk = " << maxk << endl;
	cout << "maxr = " << r << endl;
	cout << "maxa = " << maxa << endl;
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
	bigint N;
	N = 7069;
	cin >> N;
	init(N);
	if(AKS(N))
		cout << "YES" << endl;
	else
		cout << "NO" << endl;
	return 0;
}
