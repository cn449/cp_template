// 傾きが単調減少, min を返す
class CHT_monotone {
	ll L, R, M;
	vector<P<P<ll>>> segment;
	ll eval(int idx, ll x) { return segment[idx].first.first * x + segment[idx].first.second; }
public:
	CHT_monotone(ll L_, ll R_) {
		L = L_; R = R_; M = INF;
	}
	// ax + b を追加
	void add(ll a, ll b) {
		assert(M > a);
		M = a;
		while (true) {
			if (segment.empty()) break;
			ll p = segment.back().second.first;
			ll e0 = eval(ssize(segment) - 1, p);
			ll e1 = a * p + b;
			if (e0 > e1) segment.pop_back();
			else {
				ll d = segment.back().first.first - a;
				ll len = (e1 - e0 + d - 1) / d;
				segment.back().second.second = segment.back().second.first + len;
				break;
			}
		}
		ll t = L;
		if (!segment.empty()) t = segment.back().second.second;
		segment.push_back({ {a,b},{t,R} });
	}
	ll get_min(ll x) {
		assert(!segment.empty());
		assert(L <= x && x <= R);
		int small = 0, large = ssize(segment);
		while (large - small > 1) {
			int mid = midpoint(small, large);
			if (x >= segment[mid].second.first) small = mid;
			else large = mid;
		}
		return eval(small, x);
	}
	void print() {
		for (auto [p0, p1] : segment) {
			auto [a, b] = p0;
			auto [l, r] = p1;
			debug(a, b, l, r);
		}
	}
};

class CHT {
	struct Line {
		mutable ll a, b;
		mutable ll p;
		bool operator<(const Line& o) const { return a > o.a; }
		bool operator<(ll x) const { return p < x; }
	};
	ll L, R;
	multiset<Line, less<>> st;
	static ll floordiv(ll x, ll y) {
		assert(y != 0);
		if (y < 0) x = -x, y = -y;
		if (x >= 0) return x / y;
		return -((-x + y - 1) / y);
	}
	bool isect(multiset<Line>::iterator x, multiset<Line>::iterator y) {
		if (y == st.end()) { x->p = INF; return false; }
		if (x->a == y->a) x->p = (x->b <= y->b ? INF : -INF);
		else x->p = floordiv(y->b - x->b, x->a - y->a);
		return x->p >= y->p;
	}
public:
	CHT(ll L_, ll R_) : L(L_), R(R_) {}
	void add(ll a, ll b) {
		auto z = st.insert({ a, b, 0 });
		auto y = z++;
		auto x = y;
		while (isect(y, z)) z = st.erase(z);
		if (x != st.begin() && isect(--x, y)) isect(x, y = st.erase(y));
		while ((y = x) != st.begin() && (--x)->p >= y->p) isect(x, st.erase(y));
	}
	ll get_min(ll x) const {
		assert(!st.empty());
		assert(L <= x && x <= R);
		auto l = *st.lower_bound(x);
		return l.a * x + l.b;
	}
};

using T = double;
class CHT_double {
	struct Line {
		mutable T a, b;
		mutable T p;
		bool operator<(const Line& o) const { return a > o.a; }
		bool operator<(T x) const { return p < x; }
	};
	T L, R;
	multiset<Line, less<>> st;
	bool isect(multiset<Line>::iterator x, multiset<Line>::iterator y) {
		if (y == st.end()) { x->p = INF; return false; }
		if (x->a == y->a) x->p = (x->b <= y->b ? INF : -INF);
		else x->p = (y->b - x->b) / (x->a - y->a);
		return x->p >= y->p;
	}
public:
	CHT_double(T L_, T R_) : L(L_), R(R_) {}
	void add(T a, T b) {
		auto z = st.insert({ a, b, 0 });
		auto y = z++;
		auto x = y;
		while (isect(y, z)) z = st.erase(z);
		if (x != st.begin() && isect(--x, y)) isect(x, y = st.erase(y));
		while ((y = x) != st.begin() && (--x)->p >= y->p) isect(x, st.erase(y));
	}
	T get_min(T x) const {
		assert(!st.empty());
		assert(L <= x && x <= R);
		auto l = *st.lower_bound(x);
		return l.a * x + l.b;
	}
};