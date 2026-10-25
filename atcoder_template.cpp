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
using namespace std;
#define all(a) (a).begin(), (a).end()
#define rep(i, n) for (ll i = 0; i < (n); i++)
#define For(i, a, b) for (ll i = (a); i < (b); i++)
#define debug(...) cerr << __LINE__ << " | ", debug_out(#__VA_ARGS__, __VA_ARGS__)
#define test_only(x, y, ...) do { if ((x) != (y)) { debug(x, y, ##__VA_ARGS__); } } while (0)
#define test_same(x, y, ...) do { if ((x) != (y)) { debug(x, y, ##__VA_ARGS__); abort(); } } while (0)
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

// https://github.com/cn449/cp_template

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

#include <atcoder/all>
using namespace atcoder;

constexpr ll mod = mod998;
using mint = static_modint<mod>;

istream& operator >>(istream& s, mint& m) { ll y; s >> y; m = y; return s; }
istream& operator >>(istream& s, vector<mint>& v) { for (auto& e : v) { ll y; s >> y; e = y; } return s; }
ostream& operator <<(ostream& s, const mint& m) { return s << m.val(); }
ostream& operator <<(ostream& s, const vector<mint>& v) { for (auto& e : v) s << e.val() << ' '; return s; }

void debug_out(const char*) { cerr << "\n"; }
template <typename T, typename... Args> void debug_out(const char* names, T value, Args... args) {
	while (*names == ' ') ++names;
	const char* comma = strchr(names, ',');
	if (!comma) { cerr << names << ":" << value << "\n"; }
	else { cerr.write(names, comma - names) << ":" << value << ", "; debug_out(comma + 1, args...); }
}

vector<mint> fac, inv, facinv;

void modcalc(int n) {
	assert(fac.empty()); assert(inv.empty()); assert(facinv.empty());
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
	if (n < 0 or k < 0 or n < k) return 0;
	return fac[n] * facinv[k] * facinv[n - k];
}

mint perm(ll n, ll k) {
	if (n < 0 or k < 0 or n < k) return 0;
	return fac[n] * facinv[n - k];
}

mint hom(ll n, ll k) {
	if (n < 0 or k < 0 or (n == 0 && k > 0)) return 0;
	if (n == 0 && k == 0) return 1;
	return fac[n + k - 1] * facinv[k] * facinv[n - 1];
}

int main() {

	[[maybe_unused]] bool DEBUG = false;

	int TEST = 1;
	// cin >> TEST;
	while (TEST--) {
	}
}