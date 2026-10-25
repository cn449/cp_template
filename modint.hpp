constexpr ll mod = mod998;
template <int mod> class modint {
public:
	int n;
	modint() : n(0) {};
	modint(ll n_) {
		n = n_ % mod;
		if (n < 0) n += mod;
	}
	modint operator -() const { return n > 0 ? mod - n : -n; }
	bool operator == (const modint &m) const { return n == m.n; }
	bool operator != (const modint &m) const { return n != m.n; }
	modint &operator += (const modint &m) { n += m.n; if (n >= mod) n -= mod; return *this; }
	modint &operator -= (const modint &m) { n -= m.n; if (n < 0) n += mod; return *this; }
	modint &operator *= (const modint &m) { n = ll(n) * m.n % mod; return *this; }
	modint &operator /= (const modint &m) { n = ll(n) * modinv(m).n % mod; return *this; }
	modint operator +(modint m) const { return modint(*this) += m; }
	modint operator -(modint m) const { return modint(*this) -= m; }
	modint operator *(modint m) const { return modint(*this) *= m; }
	modint operator /(modint m) const { return modint(*this) /= m; }
	modint &operator ++ () { *this += 1; return *this; }
	modint operator ++ (int) { *this += 1; return *this - 1; }
	modint &operator -- () { *this -= 1; return *this; }
	modint operator -- (int) { *this -= 1; return *this + 1; }
	modint pow(ll b) const {
		modint res = 1, a = modint(*this);
		while (b) {
			if (b & 1) res *= a;
			a *= a;
			b >>= 1;
		}
		return res;
	}
	int val() const { return n; }
	friend istream &operator >> (istream &s, modint<mod> &a) { ll x; if (s >> x) a = modint(x); return s; }
	friend ostream &operator << (ostream &s, modint<mod> &a) { s << a.n; return s; }
};

using mint = modint<mod>;

mint modinv(mint x) {
	ll a = x.n;
	if (a == 0) abort();
	ll b = mod, u = 1, v = 0;
	while (b) {
		ll t = a / b;
		a -= t * b; swap(a, b);
		u -= t * v; swap(u, v);
	}
	mint res = u;
	return res;
}

vector<mint> calc_many_inverse(vector<mint> &f) {
	int n = ssize(f);
	mint prod = 1;
	rep(i, n) prod *= f[i];
	mint inv_prod = mint(1) / prod;
	vector<mint> l(n + 1), r(n + 1);
	l[0] = 1; r[n] = 1;
	rep(i, n) l[i + 1] = l[i] * f[i];
	for (int i = n - 1; i >= 0; i--) r[i] = r[i + 1] * f[i];
	vector<mint> res(n);
	rep(i, n) res[i] = inv_prod * l[i] * r[i + 1];
	return res;
}

// 0 には 0 を返す
vector<mint> calc_many_inverse_maybe_zero(vector<mint> &f) {
	int n = ssize(f);
	vector<int> idx; vector<mint> g;
	rep(i, n) {
		if (f[i] != 0) {
			idx.push_back(i);
			g.push_back(f[i]);
		}
	}
	g = calc_many_inverse(g);
	vector<mint> res(n);
	int sz = ssize(g);
	rep(i, sz) res[idx[i]] = g[i];
	return res;
}