constexpr ull mod61 = (1ULL << 61) - 1;
class mint61 {
	static const ull mask30 = (1ULL << 30) - 1;
	static const ull mask31 = (1ULL << 31) - 1;
	ull modulo(ull a) {
		ull au = a >> 61, ad = a & mod61;
		ull res = au + ad;
		if (res >= mod61) res -= mod61;
		return res;
	}
	ull multiply(ull a, ull b) {
		ull au = a >> 31, ad = a & mask31;
		ull bu = b >> 31, bd = b & mask31;
		ull mid = ad * bu + au * bd;
		ull midu = mid >> 30, midd = mid & mask30;
		return modulo(au * bu * 2 + midu + (midd << 31) + ad * bd);
	}
	ull n;
public:
	mint61() : n(0) {};
	template<integral T> mint61(T n_) {
		if constexpr (is_signed_v<T>) {
			if (n_ < 0) {
				ull tmp = ull(-(n_ + 1)) + 1;
				n = modulo(tmp);
				if (n) n = mod61 - n;
				return;
			}
		}
		n = modulo(ull(n_));
	}
	mint61 operator -() const { return mod61 - n; }
	bool operator == (const mint61& m) const { return n == m.n; }
	bool operator != (const mint61& m) const { return n != m.n; }
	mint61& operator += (const mint61& m) { n += m.n; if (n >= mod61) n -= mod61; return *this; }
	mint61& operator -= (const mint61& m) { n += mod61 - m.n; if (n >= mod61) n -= mod61; return *this; }
	mint61& operator *= (const mint61& m) { n = multiply(n, m.n); return *this; }
	mint61 operator +(mint61 m) const { return mint61(*this) += m; }
	mint61 operator -(mint61 m) const { return mint61(*this) -= m; }
	mint61 operator *(mint61 m) const { return mint61(*this) *= m; }
	mint61& operator ++ () { *this += 1; return *this; }
	mint61 operator ++ (int) { *this += 1; return *this - 1; }
	mint61& operator -- () { *this -= 1; return *this; }
	mint61 operator -- (int) { *this -= 1; return *this + 1; }
	ull val() const { return n; }
	friend istream& operator >> (istream& s, mint61& a) {
		ll x;
		if (s >> x) a = mint61(x);
		return s;
	}
	friend ostream& operator << (ostream& s, const mint61& a) { return s << a.n; }
};

template<class T> class rollinghash {
public:
	inline static const mint61 base = rng(1, mod61 - 1);
	vector<mint61> hashvalue;
	inline static vector<mint61> basepower;
	int n;
	rollinghash(const T& s) {
		n = ssize(s);
		hashvalue.resize(n + 1);
		rep(i, n) hashvalue[i + 1] = hashvalue[i] * base + (ll(s[i]) + (1ULL << 32));
		int pre = ssize(basepower);
		if (pre < n + 1) {
			basepower.resize(n + 1);
			if (pre == 0) basepower[0] = 1;
			For(i, max(0, pre - 1), n) basepower[i + 1] = basepower[i] * base;
		}
	}
	mint61 hash(ll l, ll r) const {
		assert(0 <= l && l <= r && r <= n);
		mint61 res = hashvalue[r] - hashvalue[l] * basepower[r - l];
		return res;
	}
};
