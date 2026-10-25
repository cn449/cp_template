#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <array>
#include <iomanip>
#include <utility>
#include <tuple>
#include <functional>
#include <bitset>
#include <cassert>
#include <complex>
#include <stdio.h>
#include <time.h>
#include <numeric>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <cstring>
#include <bit>
#include <chrono>
#define all(a) (a).begin(), (a).end()
#define rep(i, n) for (ll i = 0; i < (n); i++)
#define For(i, a, b) for (ll i = (a); i < (b); i++)
#define debug(...) cerr << __LINE__ << " | ", debug_out(#__VA_ARGS__, __VA_ARGS__)
#define test_only(x, y, ...) do { if ((x) != (y)) { debug(x, y, ##__VA_ARGS__); } } while (0)
#define test_same(x, y, ...) do { if ((x) != (y)) { debug(x, y, ##__VA_ARGS__); abort(); } } while (0)
using namespace std;
typedef long long ll;
typedef unsigned int uint;
typedef unsigned long long ull;
typedef long double ld;
template<class T> using P = pair<T, T>;
template<class T> using pri_l = priority_queue<T>;
template<class T> using pri_s = priority_queue<T, vector<T>, greater<T>>;
constexpr int inf = 1000000010;
constexpr ll INF = 1000000000000000010;
constexpr int mod1e9 = 1000000007;
constexpr int mod998 = 998244353;
constexpr ld eps = 1e-12;
constexpr ld pi = 3.141592653589793238;
constexpr ll ten(int n) { return n ? 10 * ten(n - 1) : 1; };
int dx[] = { 1,0,-1,0,1,1,-1,-1,0 }; int dy[] = { 0,1,0,-1,1,-1,1,-1,0 };
ll safe_mul(ll a, ll b, ll LIM = INF) { return (b != 0 && a > LIM / b ? LIM : a * b); }
void fail() { cout << "-1\n"; exit(0); } void no() { cout << "No\n"; exit(0); }
template<class T> void er(T a) { cout << a << '\n'; exit(0); }
template<class T, class U> inline bool chmax(T& a, const U& b) { if (a < b) { a = b; return true; } return false; }
template<class T, class U> inline bool chmin(T& a, const U& b) { if (a > b) { a = b; return true; } return false; }
template<class T, class U> ostream& operator << (ostream& s, const pair<T, U>& p) { s << p.first << ' ' << p.second; return s; }
template<class T> istream& operator >>(istream& s, vector<T>& v) { for (auto& e : v) s >> e; return s; }
template<class T> ostream& operator <<(ostream& s, const vector<T>& v) { for (auto& e : v) s << e << ' '; return s; }

struct fastio {
	fastio() {
		cin.tie(0); cout.tie(0);
		ios::sync_with_stdio(false);
		cout << fixed << setprecision(20);
		cerr << fixed << setprecision(20);
	}
}fastio_;

namespace rdv {
	random_device seed_gen;
	mt19937_64 engine(seed_gen());
	ll rnum(ll r) { return engine() % r; } // [0, r)
	ll rnum(ll l, ll r) { return rnum(r - l) + l; } // [l, r)
	ll rng(ll l, ll r) { return rnum(l, r + 1); } // [l, r]
	double rng01() { return engine() * pow(2, -64); }
	template<class T> void shuf(vector<T>& v) { shuffle(all(v), engine); }
	void shuf(string& s) { shuffle(all(s), engine); }
}

using namespace rdv;

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
	bool operator == (const modint& m) const { return n == m.n; }
	bool operator != (const modint& m) const { return n != m.n; }
	modint& operator += (const modint& m) { n += m.n; if (n >= mod) n -= mod; return *this; }
	modint& operator -= (const modint& m) { n -= m.n; if (n < 0) n += mod; return *this; }
	modint& operator *= (const modint& m) { n = ll(n) * m.n % mod; return *this; }
	modint& operator /= (const modint& m) { n = ll(n) * modinv(m).n % mod; return *this; }
	modint operator +(modint m) const { return modint(*this) += m; }
	modint operator -(modint m) const { return modint(*this) -= m; }
	modint operator *(modint m) const { return modint(*this) *= m; }
	modint operator /(modint m) const { return modint(*this) /= m; }
	modint& operator ++ () { *this += 1; return *this; }
	modint operator ++ (int) { *this += 1; return *this - 1; }
	modint& operator -- () { *this -= 1; return *this; }
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
	friend istream& operator >> (istream& s, modint<mod>& a) { ll x; if (s >> x) a = modint(x); return s; }
	friend ostream& operator << (ostream& s, const modint<mod>& a) { s << a.n; return s; }
};

using mint = modint<mod>;

ostream& operator <<(ostream& s, const vector<mint>& v) { for (auto& e : v) s << e.n << ' '; return s; }

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

template<class T> vector<int> compress(vector<T> &v) {
	int n = ssize(v);
	vector<T> tmp = v;
	sort(tmp.begin(), tmp.end());
	tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
	vector<int> res(n);
	for (int i = 0; i < n; i++) res[i] = lower_bound(tmp.begin(), tmp.end(), v[i]) - tmp.begin();
	return res;
}

#ifdef _MSC_VER
using lint = ll;
#else 
using lint = __int128_t;
#endif

void debug_out(const char*) { cerr << "\n"; }
template <typename T, typename... Args> void debug_out(const char* names, T value, Args... args) {
	while (*names == ' ') ++names;
	const char* comma = strchr(names, ',');
	if (!comma) { cerr << names << ":" << value << "\n"; }
	else { cerr.write(names, comma - names) << ":" << value << ", "; debug_out(comma + 1, args...); }
}

vector<mint> fac, inv, facinv;

void modcalc(int n) {
	fac.resize(n); inv.resize(n); facinv.resize(n);
	fac[0] = 1; fac[1] = 1; inv[1] = 1;
	facinv[0] = 1; facinv[1] = 1;
	for (ll i = 2; i < n; i++) {
		fac[i] = fac[i - 1] * i;
		inv[i] = -inv[mod % i] * (mod / i);
		facinv[i] = facinv[i - 1] * inv[i];
	}
}

mint comb(ll n, ll k) {
	if (n < 0 || k < 0 || n < k) return 0;
	return fac[n] * facinv[k] * facinv[n - k];
}

mint perm(ll n, ll k) {
	if (n < 0 || k < 0 || n < k) return 0;
	return fac[n] * facinv[n - k];
}

mint hom(ll n, ll k) {
	if (n < 0 || k < 0 || (n == 0 && k > 0)) return 0;
	if (n == 0 && k == 0) return 1;
	return fac[n + k - 1] * facinv[k] * facinv[n - 1];
}

// AtCoder Library
namespace atcoder {
	template <class S, auto op, auto e, class F, auto mapping, auto composition, auto id> struct lazy_segtree {
		static_assert(is_convertible_v<decltype(op), function<S(S, S)>>, "op must work as S(S, S)");
		static_assert(is_convertible_v<decltype(e), function<S()>>, "e must work as S()");
		static_assert(is_convertible_v<decltype(mapping), function<S(F, S)>>, "mapping must work as F(F, S)");
		static_assert(is_convertible_v<decltype(composition), function<F(F, F)>>, "compostiion must work as F(F, F)");
		static_assert(is_convertible_v<decltype(id), function<F()>>, "id must work as F()");
	public:
		lazy_segtree() : lazy_segtree(0) {}
		explicit lazy_segtree(int n) : lazy_segtree(vector<S>(n, e())) {}
		explicit lazy_segtree(const vector<S>& v) : _n(ssize(v)) {
			size = bit_ceil((unsigned int)(_n));
			log = countr_zero((unsigned int)size);
			d = vector<S>(2 * size, e()), lz = vector<F>(size, id());
			for (int i = 0; i < _n; i++) d[size + i] = v[i];
			for (int i = size - 1; i >= 1; i--) update(i);
		}
		void set(int p, S x) {
			assert(0 <= p && p < _n);
			p += size;
			for (int i = log; i >= 1; i--) push(p >> i);
			d[p] = x;
			for (int i = 1; i <= log; i++) update(p >> i);
		}
		S get(int p) {
			assert(0 <= p && p < _n);
			p += size;
			for (int i = log; i >= 1; i--) push(p >> i);
			return d[p];
		}
		S prod(int l, int r) {
			assert(0 <= l && l <= r && r <= _n);
			if (l == r) return e();
			l += size; r += size;
			for (int i = log; i >= 1; i--) {
				if (((l >> i) << i) != l) push(l >> i);
				if (((r >> i) << i) != r) push((r - 1) >> i);
			}
			S sml = e(), smr = e();
			while (l < r) {
				if (l & 1) sml = op(sml, d[l++]);
				if (r & 1) smr = op(d[--r], smr);
				l >>= 1; r >>= 1;
			}
			return op(sml, smr);
		}
		S all_prod() { return d[1]; }
		void apply(int p, F f) {
			assert(0 <= p && p < _n);
			p += size;
			for (int i = log; i >= 1; i--) push(p >> i);
			d[p] = mapping(f, d[p]);
			for (int i = 1; i <= log; i++) update(p >> i);
		}
		void apply(int l, int r, F f) {
			assert(0 <= l && l <= r && r <= _n);
			if (l == r) return;
			l += size; r += size;
			for (int i = log; i >= 1; i--) {
				if (((l >> i) << i) != l) push(l >> i);
				if (((r >> i) << i) != r) push((r - 1) >> i);
			}
			{
				int l2 = l, r2 = r;
				while (l < r) {
					if (l & 1) all_apply(l++, f);
					if (r & 1) all_apply(--r, f);
					l >>= 1; r >>= 1;
				}
				l = l2; r = r2;
			}
			for (int i = 1; i <= log; i++) {
				if (((l >> i) << i) != l) update(l >> i);
				if (((r >> i) << i) != r) update((r - 1) >> i);
			}
		}
		template <bool (*g)(S)> int max_right(int l) { return max_right(l, [](S x) { return g(x); }); }
		template <class G> int max_right(int l, G g) {
			assert(0 <= l && l <= _n); assert(g(e()));
			if (l == _n) return _n;
			l += size;
			for (int i = log; i >= 1; i--) push(l >> i);
			S sm = e();
			do {
				while (l % 2 == 0) l >>= 1;
				if (!g(op(sm, d[l]))) {
					while (l < size) {
						push(l);  l = 2 * l;
						if (g(op(sm, d[l]))) { sm = op(sm, d[l]); l++; }
					}
					return l - size;
				}
				sm = op(sm, d[l]); l++;
			} while ((l & -l) != l);
			return _n;
		}
		template <bool (*g)(S)> int min_left(int r) { return min_left(r, [](S x) { return g(x); }); }
		template <class G> int min_left(int r, G g) {
			assert(0 <= r && r <= _n); assert(g(e()));
			if (r == 0) return 0;
			r += size;
			for (int i = log; i >= 1; i--) push((r - 1) >> i);
			S sm = e();
			do {
				r--;
				while (r > 1 && (r % 2)) r >>= 1;
				if (!g(op(d[r], sm))) {
					while (r < size) {
						push(r); r = 2 * r + 1;
						if (g(op(d[r], sm))) { sm = op(d[r], sm); r--; }
					}
					return r + 1 - size;
				}
				sm = op(d[r], sm);
			} while ((r & -r) != r);
			return 0;
		}
	private:
		int _n, size, log;
		vector<S> d; vector<F> lz;
		void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
		void all_apply(int k, F f) {
			d[k] = mapping(f, d[k]);
			if (k < size) lz[k] = composition(f, lz[k]);
		}
		void push(int k) {
			all_apply(2 * k, lz[k]);
			all_apply(2 * k + 1, lz[k]);
			lz[k] = id();
		}
	};

	template <class S, auto op, auto e> struct segtree {
		static_assert(is_convertible_v<decltype(op), function<S(S, S)>>, "op must work as S(S, S)");
		static_assert(is_convertible_v<decltype(e), function<S()>>, "e must work as S()");
	public:
		segtree() : segtree(0) {}
		explicit segtree(int n) : segtree(vector<S>(n, e())) {}
		explicit segtree(const vector<S>& v) : _n(ssize(v)) {
			size = bit_ceil((unsigned int)(_n));
			log = countr_zero((unsigned int)size);
			d = vector<S>(2 * size, e());
			for (int i = 0; i < _n; i++) d[size + i] = v[i];
			for (int i = size - 1; i >= 1; i--) update(i);
		}
		void set(int p, S x) {
			assert(0 <= p && p < _n);
			p += size;
			d[p] = x;
			for (int i = 1; i <= log; i++) update(p >> i);
		}
		S get(int p) const { assert(0 <= p && p < _n); return d[p + size]; }
		S prod(int l, int r) const {
			assert(0 <= l && l <= r && r <= _n);
			S sml = e(), smr = e();
			l += size; r += size;
			while (l < r) {
				if (l & 1) sml = op(sml, d[l++]);
				if (r & 1) smr = op(d[--r], smr);
				l >>= 1; r >>= 1;
			}
			return op(sml, smr);
		}
		S all_prod() const { return d[1]; }
		template <bool (*f)(S)> int max_right(int l) const { return max_right(l, [](S x) { return f(x); }); }
		template <class F> int max_right(int l, F f) const {
			assert(0 <= l && l <= _n); assert(f(e()));
			if (l == _n) return _n;
			l += size;
			S sm = e();
			do {
				while (l % 2 == 0) l >>= 1;
				if (!f(op(sm, d[l]))) {
					while (l < size) {
						l = 2 * l;
						if (f(op(sm, d[l]))) { sm = op(sm, d[l]); l++; }
					}
					return l - size;
				}
				sm = op(sm, d[l]); l++;
			} while ((l & -l) != l);
			return _n;
		}
		template <bool (*f)(S)> int min_left(int r) const { return min_left(r, [](S x) { return f(x); }); }
		template <class F> int min_left(int r, F f) const {
			assert(0 <= r && r <= _n); assert(f(e()));
			if (r == 0) return 0;
			r += size;
			S sm = e();
			do {
				r--;
				while (r > 1 && (r % 2)) r >>= 1;
				if (!f(op(d[r], sm))) {
					while (r < size) {
						r = 2 * r + 1;
						if (f(op(d[r], sm))) { sm = op(d[r], sm); r--; }
					}
					return r + 1 - size;
				}
				sm = op(d[r], sm);
			} while ((r & -r) != r);
			return 0;
		}
	private:
		int _n, size, log;
		vector<S> d;
		void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
	};
}

using namespace atcoder;

int main() {

	[[maybe_unused]] bool DEBUG = false;

	int TEST = 1;
	cin >> TEST;
	while (TEST--) {
	}
}
