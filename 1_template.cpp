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

template<class T> vector<int> compress(vector<T> v) {
	int n = ssize(v);
	vector<T> tmp = v;
	sort(tmp.begin(), tmp.end());
	tmp.erase(unique(tmp.begin(), tmp.end()), tmp.end());
	vector<int> res(n);
	for (int i = 0; i < n; i++) res[i] = lower_bound(tmp.begin(), tmp.end(), v[i]) - tmp.begin();
	return res;
}

void debug_out(const char*) { cerr << "\n"; }
template <typename T, typename... Args> void debug_out(const char* names, T value, Args... args) {
	while (*names == ' ') ++names;
	const char* comma = strchr(names, ',');
	if (!comma) { cerr << names << ":" << value << "\n"; }
	else { cerr.write(names, comma - names) << ":" << value << ", "; debug_out(comma + 1, args...); }
}

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

ll modinv(ll a, ll m) {
	if (m == 1) return a;
	a %= m;
	if (a == 0) abort();
	ll b = m, u = 1, v = 0;
	while (b) {
		ll t = a / b;
		a -= t * b; swap(a, b);
		u -= t * v; swap(u, v);
	}
	u %= m;
	if (u < 0) u += m;
	return u;
}

#ifndef ATCODER_INTERNAL_MATH_HPP
#define ATCODER_INTERNAL_MATH_HPP 1

#include <utility>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace atcoder {

	namespace internal {

		// @param m `1 <= m`
		// @return x mod m
		constexpr long long safe_mod(long long x, long long m) {
			x %= m;
			if (x < 0) x += m;
			return x;
		}

		// Fast modular multiplication by barrett reduction
		// Reference: https://en.wikipedia.org/wiki/Barrett_reduction
		// NOTE: reconsider after Ice Lake
		struct barrett {
			unsigned int _m;
			unsigned long long im;

			// @param m `1 <= m`
			explicit barrett(unsigned int m) : _m(m), im((unsigned long long)(-1) / m + 1) {}

			// @return m
			unsigned int umod() const { return _m; }

			// @param a `0 <= a < m`
			// @param b `0 <= b < m`
			// @return `a * b % m`
			unsigned int mul(unsigned int a, unsigned int b) const {
				// [1] m = 1
				// a = b = im = 0, so okay

				// [2] m >= 2
				// im = ceil(2^64 / m)
				// -> im * m = 2^64 + r (0 <= r < m)
				// let z = a*b = c*m + d (0 <= c, d < m)
				// a*b * im = (c*m + d) * im = c*(im*m) + d*im = c*2^64 + c*r + d*im
				// c*r + d*im < m * m + m * im < m * m + 2^64 + m <= 2^64 + m * (m + 1) < 2^64 * 2
				// ((ab * im) >> 64) == c or c + 1
				unsigned long long z = a;
				z *= b;
#ifdef _MSC_VER
				unsigned long long x;
				_umul128(z, im, &x);
#else
				unsigned long long x =
					(unsigned long long)(((unsigned __int128)(z)*im) >> 64);
#endif
				unsigned long long y = x * _m;
				return (unsigned int)(z - y + (z < y ? _m : 0));
			}
		};

		// @param n `0 <= n`
		// @param m `1 <= m`
		// @return `(x ** n) % m`
		constexpr long long pow_mod_constexpr(long long x, long long n, int m) {
			if (m == 1) return 0;
			unsigned int _m = (unsigned int)(m);
			unsigned long long r = 1;
			unsigned long long y = safe_mod(x, m);
			while (n) {
				if (n & 1) r = (r * y) % _m;
				y = (y * y) % _m;
				n >>= 1;
			}
			return r;
		}

		// Reference:
		// M. Forisek and J. Jancina,
		// Fast Primality Testing for Integers That Fit into a Machine Word
		// @param n `0 <= n`
		constexpr bool is_prime_constexpr(int n) {
			if (n <= 1) return false;
			if (n == 2 || n == 7 || n == 61) return true;
			if (n % 2 == 0) return false;
			long long d = n - 1;
			while (d % 2 == 0) d /= 2;
			constexpr long long bases[3] = { 2, 7, 61 };
			for (long long a : bases) {
				long long t = d;
				long long y = pow_mod_constexpr(a, t, n);
				while (t != n - 1 && y != 1 && y != n - 1) {
					y = y * y % n;
					t <<= 1;
				}
				if (y != n - 1 && t % 2 == 0) {
					return false;
				}
			}
			return true;
		}
		template <int n> constexpr bool is_prime = is_prime_constexpr(n);

		// @param b `1 <= b`
		// @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b/g
		constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
			a = safe_mod(a, b);
			if (a == 0) return { b, 0 };

			// Contracts:
			// [1] s - m0 * a = 0 (mod b)
			// [2] t - m1 * a = 0 (mod b)
			// [3] s * |m1| + t * |m0| <= b
			long long s = b, t = a;
			long long m0 = 0, m1 = 1;

			while (t) {
				long long u = s / t;
				s -= t * u;
				m0 -= m1 * u;  // |m1 * u| <= |m1| * s <= b

				// [3]:
				// (s - t * u) * |m1| + t * |m0 - m1 * u|
				// <= s * |m1| - t * u * |m1| + t * (|m0| + |m1| * u)
				// = s * |m1| + t * |m0| <= b

				auto tmp = s;
				s = t;
				t = tmp;
				tmp = m0;
				m0 = m1;
				m1 = tmp;
			}
			// by [3]: |m0| <= b/g
			// by g != b: |m0| < b/g
			if (m0 < 0) m0 += b / s;
			return { s, m0 };
		}

		// Compile time primitive root
		// @param m must be prime
		// @return primitive root (and minimum in now)
		constexpr int primitive_root_constexpr(int m) {
			if (m == 2) return 1;
			if (m == 167772161) return 3;
			if (m == 469762049) return 3;
			if (m == 754974721) return 11;
			if (m == 998244353) return 3;
			int divs[20] = {};
			divs[0] = 2;
			int cnt = 1;
			int x = (m - 1) / 2;
			while (x % 2 == 0) x /= 2;
			for (int i = 3; (long long)(i)*i <= x; i += 2) {
				if (x % i == 0) {
					divs[cnt++] = i;
					while (x % i == 0) {
						x /= i;
					}
				}
			}
			if (x > 1) {
				divs[cnt++] = x;
			}
			for (int g = 2;; g++) {
				bool ok = true;
				for (int i = 0; i < cnt; i++) {
					if (pow_mod_constexpr(g, (m - 1) / divs[i], m) == 1) {
						ok = false;
						break;
					}
				}
				if (ok) return g;
			}
		}
		template <int m> constexpr int primitive_root = primitive_root_constexpr(m);

		// @param n `n < 2^32`
		// @param m `1 <= m < 2^32`
		// @return sum_{i=0}^{n-1} floor((ai + b) / m) (mod 2^64)
		unsigned long long floor_sum_unsigned(unsigned long long n,
			unsigned long long m,
			unsigned long long a,
			unsigned long long b) {
			unsigned long long ans = 0;
			while (true) {
				if (a >= m) {
					ans += n * (n - 1) / 2 * (a / m);
					a %= m;
				}
				if (b >= m) {
					ans += n * (b / m);
					b %= m;
				}

				unsigned long long y_max = a * n + b;
				if (y_max < m) break;
				// y_max < m * (n + 1)
				// floor(y_max / m) <= n
				n = (unsigned long long)(y_max / m);
				b = (unsigned long long)(y_max % m);
				std::swap(m, a);
			}
			return ans;
		}

	}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_INTERNAL_MATH_HPP

#ifndef ATCODER_INTERNAL_TYPE_TRAITS_HPP
#define ATCODER_INTERNAL_TYPE_TRAITS_HPP 1

#include <cassert>
#include <numeric>
#include <type_traits>

namespace atcoder {

	namespace internal {

#ifndef _MSC_VER
		template <class T>
		using is_signed_int128 =
			typename std::conditional<std::is_same<T, __int128_t>::value ||
			std::is_same<T, __int128>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using is_unsigned_int128 =
			typename std::conditional<std::is_same<T, __uint128_t>::value ||
			std::is_same<T, unsigned __int128>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using make_unsigned_int128 =
			typename std::conditional<std::is_same<T, __int128_t>::value,
			__uint128_t,
			unsigned __int128>;

		template <class T>
		using is_integral = typename std::conditional<std::is_integral<T>::value ||
			is_signed_int128<T>::value ||
			is_unsigned_int128<T>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using is_signed_int = typename std::conditional<(is_integral<T>::value&&
			std::is_signed<T>::value) ||
			is_signed_int128<T>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using is_unsigned_int =
			typename std::conditional<(is_integral<T>::value&&
				std::is_unsigned<T>::value) ||
			is_unsigned_int128<T>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using to_unsigned = typename std::conditional<
			is_signed_int128<T>::value,
			make_unsigned_int128<T>,
			typename std::conditional<std::is_signed<T>::value,
			std::make_unsigned<T>,
			std::common_type<T>>::type>::type;

#else

		template <class T> using is_integral = typename std::is_integral<T>;

		template <class T>
		using is_signed_int =
			typename std::conditional<is_integral<T>::value&& std::is_signed<T>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using is_unsigned_int =
			typename std::conditional<is_integral<T>::value&&
			std::is_unsigned<T>::value,
			std::true_type,
			std::false_type>::type;

		template <class T>
		using to_unsigned = typename std::conditional<is_signed_int<T>::value,
			std::make_unsigned<T>,
			std::common_type<T>>::type;

#endif

		template <class T>
		using is_signed_int_t = std::enable_if_t<is_signed_int<T>::value>;

		template <class T>
		using is_unsigned_int_t = std::enable_if_t<is_unsigned_int<T>::value>;

		template <class T> using to_unsigned_t = typename to_unsigned<T>::type;

	}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_INTERNAL_TYPE_TRAITS_HPP

#ifndef ATCODER_MODINT_HPP
#define ATCODER_MODINT_HPP 1

#include <cassert>
#include <numeric>
#include <type_traits>

#ifdef _MSC_VER
#include <intrin.h>
#endif

// #include "atcoder/internal_math"
// #include "atcoder/internal_type_traits"

namespace atcoder {

	namespace internal {

		struct modint_base {};
		struct static_modint_base : modint_base {};

		template <class T> using is_modint = std::is_base_of<modint_base, T>;
		template <class T> using is_modint_t = std::enable_if_t<is_modint<T>::value>;

	}  // namespace internal

	template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
	struct static_modint : internal::static_modint_base {
		using mint = static_modint;

	public:
		static constexpr int mod() { return m; }
		static mint raw(int v) {
			mint x;
			x._v = v;
			return x;
		}

		static_modint() : _v(0) {}
		template <class T, internal::is_signed_int_t<T>* = nullptr>
		static_modint(T v) {
			long long x = (long long)(v % (long long)(umod()));
			if (x < 0) x += umod();
			_v = (unsigned int)(x);
		}
		template <class T, internal::is_unsigned_int_t<T>* = nullptr>
		static_modint(T v) {
			_v = (unsigned int)(v % umod());
		}

		int val() const { return _v; }

		mint& operator++() {
			_v++;
			if (_v == umod()) _v = 0;
			return *this;
		}
		mint& operator--() {
			if (_v == 0) _v = umod();
			_v--;
			return *this;
		}
		mint operator++(int) {
			mint result = *this;
			++* this;
			return result;
		}
		mint operator--(int) {
			mint result = *this;
			--* this;
			return result;
		}

		mint& operator+=(const mint& rhs) {
			_v += rhs._v;
			if (_v >= umod()) _v -= umod();
			return *this;
		}
		mint& operator-=(const mint& rhs) {
			_v -= rhs._v;
			if (_v >= umod()) _v += umod();
			return *this;
		}
		mint& operator*=(const mint& rhs) {
			unsigned long long z = _v;
			z *= rhs._v;
			_v = (unsigned int)(z % umod());
			return *this;
		}
		mint& operator/=(const mint& rhs) { return *this = *this * rhs.inv(); }

		mint operator+() const { return *this; }
		mint operator-() const { return mint() - *this; }

		mint pow(long long n) const {
			assert(0 <= n);
			mint x = *this, r = 1;
			while (n) {
				if (n & 1) r *= x;
				x *= x;
				n >>= 1;
			}
			return r;
		}
		mint inv() const {
			if (prime) {
				assert(_v);
				return pow(umod() - 2);
			}
			else {
				auto eg = internal::inv_gcd(_v, m);
				assert(eg.first == 1);
				return eg.second;
			}
		}

		friend mint operator+(const mint& lhs, const mint& rhs) {
			return mint(lhs) += rhs;
		}
		friend mint operator-(const mint& lhs, const mint& rhs) {
			return mint(lhs) -= rhs;
		}
		friend mint operator*(const mint& lhs, const mint& rhs) {
			return mint(lhs) *= rhs;
		}
		friend mint operator/(const mint& lhs, const mint& rhs) {
			return mint(lhs) /= rhs;
		}
		friend bool operator==(const mint& lhs, const mint& rhs) {
			return lhs._v == rhs._v;
		}
		friend bool operator!=(const mint& lhs, const mint& rhs) {
			return lhs._v != rhs._v;
		}

	private:
		unsigned int _v;
		static constexpr unsigned int umod() { return m; }
		static constexpr bool prime = internal::is_prime<m>;
	};

	template <int id> struct dynamic_modint : internal::modint_base {
		using mint = dynamic_modint;

	public:
		static int mod() { return (int)(bt.umod()); }
		static void set_mod(int m) {
			assert(1 <= m);
			bt = internal::barrett(m);
		}
		static mint raw(int v) {
			mint x;
			x._v = v;
			return x;
		}

		dynamic_modint() : _v(0) {}
		template <class T, internal::is_signed_int_t<T>* = nullptr>
		dynamic_modint(T v) {
			long long x = (long long)(v % (long long)(mod()));
			if (x < 0) x += mod();
			_v = (unsigned int)(x);
		}
		template <class T, internal::is_unsigned_int_t<T>* = nullptr>
		dynamic_modint(T v) {
			_v = (unsigned int)(v % mod());
		}

		int val() const { return _v; }

		mint& operator++() {
			_v++;
			if (_v == umod()) _v = 0;
			return *this;
		}
		mint& operator--() {
			if (_v == 0) _v = umod();
			_v--;
			return *this;
		}
		mint operator++(int) {
			mint result = *this;
			++* this;
			return result;
		}
		mint operator--(int) {
			mint result = *this;
			--* this;
			return result;
		}

		mint& operator+=(const mint& rhs) {
			_v += rhs._v;
			if (_v >= umod()) _v -= umod();
			return *this;
		}
		mint& operator-=(const mint& rhs) {
			_v += mod() - rhs._v;
			if (_v >= umod()) _v -= umod();
			return *this;
		}
		mint& operator*=(const mint& rhs) {
			_v = bt.mul(_v, rhs._v);
			return *this;
		}
		mint& operator/=(const mint& rhs) { return *this = *this * rhs.inv(); }

		mint operator+() const { return *this; }
		mint operator-() const { return mint() - *this; }

		mint pow(long long n) const {
			assert(0 <= n);
			mint x = *this, r = 1;
			while (n) {
				if (n & 1) r *= x;
				x *= x;
				n >>= 1;
			}
			return r;
		}
		mint inv() const {
			auto eg = internal::inv_gcd(_v, mod());
			assert(eg.first == 1);
			return eg.second;
		}

		friend mint operator+(const mint& lhs, const mint& rhs) {
			return mint(lhs) += rhs;
		}
		friend mint operator-(const mint& lhs, const mint& rhs) {
			return mint(lhs) -= rhs;
		}
		friend mint operator*(const mint& lhs, const mint& rhs) {
			return mint(lhs) *= rhs;
		}
		friend mint operator/(const mint& lhs, const mint& rhs) {
			return mint(lhs) /= rhs;
		}
		friend bool operator==(const mint& lhs, const mint& rhs) {
			return lhs._v == rhs._v;
		}
		friend bool operator!=(const mint& lhs, const mint& rhs) {
			return lhs._v != rhs._v;
		}

	private:
		unsigned int _v;
		static internal::barrett bt;
		static unsigned int umod() { return bt.umod(); }
	};
	template <int id> internal::barrett dynamic_modint<id>::bt(998244353);

	using modint998244353 = static_modint<998244353>;
	using modint1000000007 = static_modint<1000000007>;
	using modint = dynamic_modint<-1>;

	namespace internal {

		template <class T>
		using is_static_modint = std::is_base_of<internal::static_modint_base, T>;

		template <class T>
		using is_static_modint_t = std::enable_if_t<is_static_modint<T>::value>;

		template <class> struct is_dynamic_modint : public std::false_type {};
		template <int id>
		struct is_dynamic_modint<dynamic_modint<id>> : public std::true_type {};

		template <class T>
		using is_dynamic_modint_t = std::enable_if_t<is_dynamic_modint<T>::value>;

	}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_MODINT_HPP

#ifndef ATCODER_INTERNAL_BITOP_HPP
#define ATCODER_INTERNAL_BITOP_HPP 1

#ifdef _MSC_VER
#include <intrin.h>
#endif

#if __cplusplus >= 202002L
#include <bit>
#endif

namespace atcoder {

	namespace internal {

#if __cplusplus >= 202002L

		using std::bit_ceil;

#else

		// @return same with std::bit::bit_ceil
		unsigned int bit_ceil(unsigned int n) {
			unsigned int x = 1;
			while (x < (unsigned int)(n)) x *= 2;
			return x;
		}

#endif

		// @param n `1 <= n`
		// @return same with std::bit::countr_zero
		int countr_zero(unsigned int n) {
#ifdef _MSC_VER
			unsigned long index;
			_BitScanForward(&index, n);
			return index;
#else
			return __builtin_ctz(n);
#endif
		}

		// @param n `1 <= n`
		// @return same with std::bit::countr_zero
		constexpr int countr_zero_constexpr(unsigned int n) {
			int x = 0;
			while (!(n & (1 << x))) x++;
			return x;
		}

	}  // namespace internal

}  // namespace atcoder

#endif  // ATCODER_INTERNAL_BITOP_HPP

#ifndef ATCODER_CONVOLUTION_HPP
#define ATCODER_CONVOLUTION_HPP 1

#include <algorithm>
#include <array>
#include <cassert>
#include <type_traits>
#include <vector>

// #include "atcoder/internal_bit"
// #include "atcoder/modint"

namespace atcoder {

	namespace internal {

		template <class mint,
			int g = internal::primitive_root<mint::mod()>,
			internal::is_static_modint_t<mint>* = nullptr>
			struct fft_info {
			static constexpr int rank2 = countr_zero_constexpr(mint::mod() - 1);
			std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
			std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1

			std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
			std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;

			std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
			std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;

			fft_info() {
				root[rank2] = mint(g).pow((mint::mod() - 1) >> rank2);
				iroot[rank2] = root[rank2].inv();
				for (int i = rank2 - 1; i >= 0; i--) {
					root[i] = root[i + 1] * root[i + 1];
					iroot[i] = iroot[i + 1] * iroot[i + 1];
				}

				{
					mint prod = 1, iprod = 1;
					for (int i = 0; i <= rank2 - 2; i++) {
						rate2[i] = root[i + 2] * prod;
						irate2[i] = iroot[i + 2] * iprod;
						prod *= iroot[i + 2];
						iprod *= root[i + 2];
					}
				}
				{
					mint prod = 1, iprod = 1;
					for (int i = 0; i <= rank2 - 3; i++) {
						rate3[i] = root[i + 3] * prod;
						irate3[i] = iroot[i + 3] * iprod;
						prod *= iroot[i + 3];
						iprod *= root[i + 3];
					}
				}
			}
		};

		template <class mint, internal::is_static_modint_t<mint>* = nullptr>
		void butterfly(std::vector<mint>& a) {
			int n = int(a.size());
			int h = internal::countr_zero((unsigned int)n);

			static const fft_info<mint> info;

			int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
			while (len < h) {
				if (h - len == 1) {
					int p = 1 << (h - len - 1);
					mint rot = 1;
					for (int s = 0; s < (1 << len); s++) {
						int offset = s << (h - len);
						for (int i = 0; i < p; i++) {
							auto l = a[i + offset];
							auto r = a[i + offset + p] * rot;
							a[i + offset] = l + r;
							a[i + offset + p] = l - r;
						}
						if (s + 1 != (1 << len))
							rot *= info.rate2[countr_zero(~(unsigned int)(s))];
					}
					len++;
				}
				else {
					// 4-base
					int p = 1 << (h - len - 2);
					mint rot = 1, imag = info.root[2];
					for (int s = 0; s < (1 << len); s++) {
						mint rot2 = rot * rot;
						mint rot3 = rot2 * rot;
						int offset = s << (h - len);
						for (int i = 0; i < p; i++) {
							auto mod2 = 1ULL * mint::mod() * mint::mod();
							auto a0 = 1ULL * a[i + offset].val();
							auto a1 = 1ULL * a[i + offset + p].val() * rot.val();
							auto a2 = 1ULL * a[i + offset + 2 * p].val() * rot2.val();
							auto a3 = 1ULL * a[i + offset + 3 * p].val() * rot3.val();
							auto a1na3imag =
								1ULL * mint(a1 + mod2 - a3).val() * imag.val();
							auto na2 = mod2 - a2;
							a[i + offset] = a0 + a2 + a1 + a3;
							a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
							a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
							a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
						}
						if (s + 1 != (1 << len))
							rot *= info.rate3[countr_zero(~(unsigned int)(s))];
					}
					len += 2;
				}
			}
		}

		template <class mint, internal::is_static_modint_t<mint>* = nullptr>
		void butterfly_inv(std::vector<mint>& a) {
			int n = int(a.size());
			int h = internal::countr_zero((unsigned int)n);

			static const fft_info<mint> info;

			int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
			while (len) {
				if (len == 1) {
					int p = 1 << (h - len);
					mint irot = 1;
					for (int s = 0; s < (1 << (len - 1)); s++) {
						int offset = s << (h - len + 1);
						for (int i = 0; i < p; i++) {
							auto l = a[i + offset];
							auto r = a[i + offset + p];
							a[i + offset] = l + r;
							a[i + offset + p] =
								(unsigned long long)((unsigned int)(l.val() - r.val()) + mint::mod()) *
								irot.val();
							;
						}
						if (s + 1 != (1 << (len - 1)))
							irot *= info.irate2[countr_zero(~(unsigned int)(s))];
					}
					len--;
				}
				else {
					// 4-base
					int p = 1 << (h - len);
					mint irot = 1, iimag = info.iroot[2];
					for (int s = 0; s < (1 << (len - 2)); s++) {
						mint irot2 = irot * irot;
						mint irot3 = irot2 * irot;
						int offset = s << (h - len + 2);
						for (int i = 0; i < p; i++) {
							auto a0 = 1ULL * a[i + offset + 0 * p].val();
							auto a1 = 1ULL * a[i + offset + 1 * p].val();
							auto a2 = 1ULL * a[i + offset + 2 * p].val();
							auto a3 = 1ULL * a[i + offset + 3 * p].val();

							auto a2na3iimag =
								1ULL *
								mint((mint::mod() + a2 - a3) * iimag.val()).val();

							a[i + offset] = a0 + a1 + a2 + a3;
							a[i + offset + 1 * p] =
								(a0 + (mint::mod() - a1) + a2na3iimag) * irot.val();
							a[i + offset + 2 * p] =
								(a0 + a1 + (mint::mod() - a2) + (mint::mod() - a3)) *
								irot2.val();
							a[i + offset + 3 * p] =
								(a0 + (mint::mod() - a1) + (mint::mod() - a2na3iimag)) *
								irot3.val();
						}
						if (s + 1 != (1 << (len - 2)))
							irot *= info.irate3[countr_zero(~(unsigned int)(s))];
					}
					len -= 2;
				}
			}
		}

		template <class mint, internal::is_static_modint_t<mint>* = nullptr>
		std::vector<mint> convolution_naive(const std::vector<mint>& a,
			const std::vector<mint>& b) {
			int n = int(a.size()), m = int(b.size());
			std::vector<mint> ans(n + m - 1);
			if (n < m) {
				for (int j = 0; j < m; j++) {
					for (int i = 0; i < n; i++) {
						ans[i + j] += a[i] * b[j];
					}
				}
			}
			else {
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++) {
						ans[i + j] += a[i] * b[j];
					}
				}
			}
			return ans;
		}

		template <class mint, internal::is_static_modint_t<mint>* = nullptr>
		std::vector<mint> convolution_fft(std::vector<mint> a, std::vector<mint> b) {
			int n = int(a.size()), m = int(b.size());
			int z = (int)internal::bit_ceil((unsigned int)(n + m - 1));
			a.resize(z);
			internal::butterfly(a);
			b.resize(z);
			internal::butterfly(b);
			for (int i = 0; i < z; i++) {
				a[i] *= b[i];
			}
			internal::butterfly_inv(a);
			a.resize(n + m - 1);
			mint iz = mint(z).inv();
			for (int i = 0; i < n + m - 1; i++) a[i] *= iz;
			return a;
		}

	}  // namespace internal

	template <class mint, internal::is_static_modint_t<mint>* = nullptr>
	std::vector<mint> convolution(std::vector<mint>&& a, std::vector<mint>&& b) {
		int n = int(a.size()), m = int(b.size());
		if (!n || !m) return {};

		int z = (int)internal::bit_ceil((unsigned int)(n + m - 1));
		assert((mint::mod() - 1) % z == 0);

		if (std::min(n, m) <= 60) return convolution_naive(std::move(a), std::move(b));
		return internal::convolution_fft(std::move(a), std::move(b));
	}
	template <class mint, internal::is_static_modint_t<mint>* = nullptr>
	std::vector<mint> convolution(const std::vector<mint>& a,
		const std::vector<mint>& b) {
		int n = int(a.size()), m = int(b.size());
		if (!n || !m) return {};

		int z = (int)internal::bit_ceil((unsigned int)(n + m - 1));
		assert((mint::mod() - 1) % z == 0);

		if (std::min(n, m) <= 60) return convolution_naive(a, b);
		return internal::convolution_fft(a, b);
	}

	template <unsigned int mod = 998244353,
		class T,
		std::enable_if_t<internal::is_integral<T>::value>* = nullptr>
		std::vector<T> convolution(const std::vector<T>& a, const std::vector<T>& b) {
		int n = int(a.size()), m = int(b.size());
		if (!n || !m) return {};

		using mint = static_modint<mod>;

		int z = (int)internal::bit_ceil((unsigned int)(n + m - 1));
		assert((mint::mod() - 1) % z == 0);

		std::vector<mint> a2(n), b2(m);
		for (int i = 0; i < n; i++) {
			a2[i] = mint(a[i]);
		}
		for (int i = 0; i < m; i++) {
			b2[i] = mint(b[i]);
		}
		auto c2 = convolution(std::move(a2), std::move(b2));
		std::vector<T> c(n + m - 1);
		for (int i = 0; i < n + m - 1; i++) {
			c[i] = c2[i].val();
		}
		return c;
	}

	std::vector<long long> convolution_ll(const std::vector<long long>& a,
		const std::vector<long long>& b) {
		int n = int(a.size()), m = int(b.size());
		if (!n || !m) return {};

		static constexpr unsigned long long MOD1 = 754974721;  // 2^24
		static constexpr unsigned long long MOD2 = 167772161;  // 2^25
		static constexpr unsigned long long MOD3 = 469762049;  // 2^26
		static constexpr unsigned long long M2M3 = MOD2 * MOD3;
		static constexpr unsigned long long M1M3 = MOD1 * MOD3;
		static constexpr unsigned long long M1M2 = MOD1 * MOD2;
		static constexpr unsigned long long M1M2M3 = MOD1 * MOD2 * MOD3;

		static constexpr unsigned long long i1 =
			internal::inv_gcd(MOD2 * MOD3, MOD1).second;
		static constexpr unsigned long long i2 =
			internal::inv_gcd(MOD1 * MOD3, MOD2).second;
		static constexpr unsigned long long i3 =
			internal::inv_gcd(MOD1 * MOD2, MOD3).second;

		static constexpr int MAX_AB_BIT = 24;
		static_assert(MOD1 % (1ull << MAX_AB_BIT) == 1, "MOD1 isn't enough to support an array length of 2^24.");
		static_assert(MOD2 % (1ull << MAX_AB_BIT) == 1, "MOD2 isn't enough to support an array length of 2^24.");
		static_assert(MOD3 % (1ull << MAX_AB_BIT) == 1, "MOD3 isn't enough to support an array length of 2^24.");
		assert(n + m - 1 <= (1 << MAX_AB_BIT));

		auto c1 = convolution<MOD1>(a, b);
		auto c2 = convolution<MOD2>(a, b);
		auto c3 = convolution<MOD3>(a, b);

		std::vector<long long> c(n + m - 1);
		for (int i = 0; i < n + m - 1; i++) {
			unsigned long long x = 0;
			x += (c1[i] * i1) % MOD1 * M2M3;
			x += (c2[i] * i2) % MOD2 * M1M3;
			x += (c3[i] * i3) % MOD3 * M1M2;
			// B = 2^63, -B <= x, r(real value) < B
			// (x, x - M, x - 2M, or x - 3M) = r (mod 2B)
			// r = c1[i] (mod MOD1)
			// focus on MOD1
			// r = x, x - M', x - 2M', x - 3M' (M' = M % 2^64) (mod 2B)
			// r = x,
			//     x - M' + (0 or 2B),
			//     x - 2M' + (0, 2B or 4B),
			//     x - 3M' + (0, 2B, 4B or 6B) (without mod!)
			// (r - x) = 0, (0)
			//           - M' + (0 or 2B), (1)
			//           -2M' + (0 or 2B or 4B), (2)
			//           -3M' + (0 or 2B or 4B or 6B) (3) (mod MOD1)
			// we checked that
			//   ((1) mod MOD1) mod 5 = 2
			//   ((2) mod MOD1) mod 5 = 3
			//   ((3) mod MOD1) mod 5 = 4
			long long diff =
				c1[i] - internal::safe_mod((long long)(x), (long long)(MOD1));
			if (diff < 0) diff += MOD1;
			static constexpr unsigned long long offset[5] = {
				0, 0, M1M2M3, 2 * M1M2M3, 3 * M1M2M3 };
			x -= offset[diff % 5];
			c[i] = x;
		}

		return c;
	}

}  // namespace atcoder

#endif  // ATCODER_CONVOLUTION_HPP

using mint = modint;
ll mod;

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

using fps = vector<mint>;
using Poly = vector<mint>;

// 愚直掛け算
fps naive_conv(fps f, fps g) {
	int n = ssize(f), m = ssize(g);
	fps res(n + m - 1);
	rep(i, n) rep(j, m) res[i + j] += f[i] * g[j];
	return res;
}

// ABC245D
// 多項式としての割り算(愚直)
// (商, 余り) の組を返す
P<Poly> naive_div_poly(Poly f, Poly g) {
	while (!g.empty() && g.back() == 0) g.pop_back();
	assert(!g.empty());
	int n = ssize(f), m = ssize(g);
	n--; m--; // deg f = n, deg g = m
	mint inv = mint(1) / mint(g.back());
	int sz = n - m;
	vector<mint> q(sz + 1);
	for (int i = n; i >= m; i--) {
		q[i - m] = f[i] * inv;
		rep(j, m + 1) f[i - m + j] -= q[i - m] * g[j];
	}
	vector<mint> r(m);
	rep(i, m) r[i] = f[i];
	return { q, r };
}

// n 次多項式, 0,1,...,n での値が分かっている
// t での値を返す
// ARC033-D
mint interpolate(int n, vector<mint> a, ll t) {
	assert(ssize(a) >= n + 1);
	vector<mint> p(n + 1);
	rep(i, n + 1) {
		p[i] = a[i] * facinv[i] * facinv[n - i];
		if ((n - i) % 2 == 1) p[i] *= -1;
	}
	mint res = 0;
	vector<mint> l(n + 2), r(n + 2);
	l[0] = r[0] = 1;
	rep(i, n + 1) l[i + 1] = l[i] * (t - i);
	rep(i, n + 1) r[i + 1] = r[i] * (t - n + i);
	reverse(all(r));
	rep(i, n + 1) res += l[i] * r[i + 1] * p[i];
	return res;
}

// n 次多項式, 0,1,...,n での値が分かっている
// 係数列を返す, NTT を使わない
// CF901-E
vector<mint> interpolate_coeff_naive(int n, vector<mint> a) {
	assert(ssize(a) >= n + 1);
	// Σ(p[i] * (x - 0)(x - 1)...(x - n) / (x - i))
	vector<mint> p(n + 1);
	rep(i, n + 1) {
		p[i] = a[i] * facinv[i] * facinv[n - i];
		if ((n - i) % 2 == 1) p[i] *= -1;
	}
	vector<mint> res(n + 1);
	vector<mint> prod = { 1 };
	rep(i, n + 1) {
		vector<mint> f = { -i,1 };
		prod = naive_conv(prod, f);
	}
	rep(i, n + 1) {
		vector<mint> f = { -i,1 };
		auto [g, _] = naive_div_poly(prod, f);
		rep(j, n + 1) res[j] += g[j] * p[i];
	}
	return res;
}

int main() {

	[[maybe_unused]] bool DEBUG = false;

	int TEST = 1;
	cin >> TEST;
	while (TEST--) {
	}
}