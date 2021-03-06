#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

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
	ll deg = P->degree;
	for(ll i=deg;i>0;i--)
		cout << P->coeff[i] << "x^" << i << " + ";
	cout << P->coeff[0] << endl;
}

void polymod(poly* P,ll r,ll n)
{
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

void polymulmod(poly* P,poly* Q,ll r,ll n)
{
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

void polypow(poly *P,ll r,ll n)
{
	ll b=n;
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

bool isequal(poly* P,poly* Q)
{
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
	lhs->degree = 1;
	lhs->coeff[1] = 1;
	lhs->coeff[0] = a;
	polypow(lhs,r,n);

	rhs->degree = 1;
	rhs->coeff[1] = 1;
	rhs->coeff[0] = 0;
	polypow(rhs,r,n);
	rhs->coeff[0] += (a%n);
	if(rhs->coeff[0] >=n )
		rhs->coeff[0] -= n;
	polymod(rhs,r,n);
	return isequal(lhs,rhs);
}

bool AKS(ll N)
{
	if(IsPerfectPower(N))
		return false;

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

	ll maxr = r;
	if(r>=N)
		maxr = N-1;
	for(a=maxr;a>1;a--)
	{
		gcdan = gcd(a,N);
		if(gcdan!=1)
			return false;
	}

	if(N<=r)
		return true;

	ll maxa = ceil((double)sqrt(phi(r))*(double)log2(N));
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
	ll n;
	cin >> n;
	if(AKS(n))
		cout << n << " is prime." << endl;
	else
		cout << n << " is not prime." << endl;
	return 0;
}
