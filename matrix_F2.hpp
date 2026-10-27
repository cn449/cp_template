class Bitset {
	int n, sz;
	vector<ull> data;
public:
	Bitset(int n_) { n = n_; sz = (n + 63) / 64; data.resize(sz); }
	void set(int pos) { data[pos >> 6] |= 1ULL << (pos & 63); }
	void reset(int pos) { data[pos >> 6] &= ~(1ULL << (pos & 63)); }
	void flip(int pos) { data[pos >> 6] ^= 1ULL << (pos & 63); }
	bool get(int pos) const { return data[pos >> 6] >> (pos & 63) & 1; }
	void set(int l, int r) {
		if (l == r) return;
		int x = l >> 6, y = (r - 1) >> 6;
		ull a = ~0ULL << (l & 63), b = ~0ULL >> (63 - ((r - 1) & 63));
		if (x == y) { data[x] |= a & b; return; }
		data[x] |= a;
		for (int i = x + 1; i < y; i++) data[i] = ~0ULL;
		data[y] |= b;
	}
	void reset(int l, int r) {
		if (l == r) return;
		int x = l >> 6, y = (r - 1) >> 6;
		ull a = ~0ULL << (l & 63), b = ~0ULL >> (63 - ((r - 1) & 63));
		if (x == y) { data[x] &= ~(a & b); return; }
		data[x] &= ~a;
		for (int i = x + 1; i < y; i++) data[i] = 0;
		data[y] &= ~b;
	}
	void flip(int l, int r) {
		if (l == r) return;
		int x = l >> 6, y = (r - 1) >> 6;
		ull a = ~0ULL << (l & 63), b = ~0ULL >> (63 - ((r - 1) & 63));
		if (x == y) { data[x] ^= a & b; return; }
		data[x] ^= a;
		for (int i = x + 1; i < y; i++) data[i] = ~data[i];
		data[y] ^= b;
	}
	string to_string() const { string res = ""; rep(i, n) res += get(i) ? '1' : '0'; return res; }
	int count() const { int res = 0; rep(i, sz) res += popcount(data[i]); return res; }
	int count(int l, int r) const {
		if (l == r) return 0;
		int x = l >> 6, y = (r - 1) >> 6;
		ull a = ~0ULL << (l & 63), b = ~0ULL >> (63 - ((r - 1) & 63));
		if (x == y) return popcount(data[x] & a & b);
		int res = popcount(data[x] & a) + popcount(data[y] & b);
		for (int i = x + 1; i < y; i++) res += popcount(data[i]);
		return res;
	}
	int size() const { return n; }
	bool operator == (const Bitset& b) const { return n == b.n && data == b.data; }
	bool operator != (const Bitset& b) const { return !(*this == b); }
	bool operator < (const Bitset& b) const {
		rep(i, sz) if (data[i] != b.data[i]) {
			ull x = data[i] ^ b.data[i];
			ull w = x & ull(~x + 1);
			return (data[i] & w) == 0;
		}
		return false;
	}
	Bitset& operator ^= (const Bitset& b) { rep(i, sz) data[i] ^= b.data[i]; return *this; }
	Bitset& operator &= (const Bitset& b) { rep(i, sz) data[i] &= b.data[i]; return *this; }
	Bitset& operator |= (const Bitset& b) { rep(i, sz) data[i] |= b.data[i]; return *this; }
	Bitset operator ^(const Bitset& b) const { Bitset r(*this); r ^= b; return r; }
	Bitset operator &(const Bitset& b) const { Bitset r(*this); r &= b; return r; }
	Bitset operator |(const Bitset& b) const { Bitset r(*this); r |= b; return r; }
	Bitset operator ~ () const { Bitset r(*this); r.flip(0, n); return r; }
	Bitset& operator <<= (int k) {
		assert(k >= 0);
		if (k == 0) return *this;
		if (k >= n) { fill(data.begin(), data.end(), 0); return *this; }
		int a = k >> 6, b = k & 63;
		if (b) {
			for (int i = sz - 1; i > a; i--) data[i] = (data[i - a] << b) | (data[i - a - 1] >> (64 - b));
			data[a] = data[0] << b;
		}
		else move_backward(data.begin(), data.end() - a, data.end());
		fill(data.begin(), data.begin() + a, 0);
		if (n & 63) data.back() &= (1ULL << (n & 63)) - 1;
		return *this;
	}
	Bitset& operator >>= (int k) {
		assert(k >= 0);
		if (k == 0) return *this;
		if (k >= n) { fill(data.begin(), data.end(), 0); return *this; }
		int a = k >> 6, b = k & 63;
		if (b) {
			for (int i = 0; i < sz - a - 1; i++) data[i] = (data[i + a] >> b) | (data[i + a + 1] << (64 - b));
			data[sz - a - 1] = data[sz - 1] >> b;
		}
		else move(data.begin() + a, data.end(), data.begin());
		fill(data.end() - a, data.end(), 0);
		return *this;
	}
	Bitset operator << (int k) const { Bitset r(*this); r <<= k; return r; }
	Bitset operator >> (int k) const { Bitset r(*this); r >>= k; return r; }
	void or_shift(int k) {
		assert(k >= 0);
		if (k == 0 || k >= n) return;
		int a = k >> 6, b = k & 63;
		if (b) {
			for (int i = sz - 1; i > a; i--) data[i] |= (data[i - a] << b) | (data[i - a - 1] >> (64 - b));
			data[a] |= data[0] << b;
		}
		else for (int i = sz - 1; i >= a; i--) data[i] |= data[i - a];
		if (n & 63) data.back() &= (1ULL << (n & 63)) - 1;
	}
	int find_next(int pos) const {
		if (pos >= n - 1) return n;
		pos++;
		int i = pos >> 6;
		ull x = data[i] & (~0ULL << (pos & 63));
		if (x) return (i << 6) + countr_zero(x);
		for (i++; i < sz; i++) if (data[i]) return (i << 6) + countr_zero(data[i]);
		return n;
	}
	int find_first() const { return find_next(-1); }
	int find_prev(int pos) const {
		if (pos <= 0) return n;
		pos--;
		int i = pos >> 6;
		ull x = data[i] & (~0ULL >> (63 - (pos & 63)));
		if (x) return (i << 6) + 63 - countl_zero(x);
		for (i--; i >= 0; i--) if (data[i]) return (i << 6) + 63 - countl_zero(data[i]);
		return n;
	}
	int find_last() const { return find_prev(n); }
	bool operator [](int pos) const { return get(pos); }
};

using BS = Bitset;

vector<BS> mat_mul_F2(vector<BS> &a, vector<BS> &b) {
	int n = ssize(a), m = a[0].size();
	int m_ = ssize(b), k = b[0].size();
	assert(m == m_);
	vector res(n, BS(k));
	rep(i, n) rep(j, m) if (a[i][j]) res[i] ^= b[j];
	return res;
}

int mat_rank_F2(vector<BS>& a) {
	vector<BS> basis;
	for (auto e : a) {
		for (auto b : basis) if ((e ^ b) < e) e ^= b;
		if (e.count() != 0) basis.push_back(e);
	}
	return ssize(basis);
}

// O(nm^2/w) 時間
// ax = v の解 x を返す
// first : 解空間の次元（存在しないとき -1）
// second : 0 番目の要素は解の 1 つ、1 番目以降は解空間の基底
// 行とか列が 0 になるときに気を付ける（別途処理が必要がち？）
template<class T> pair<int, vector<BS>> linear_equation_F2(vector<BS> a, vector<T> v, bool unique = false) {
	int n = ssize(a), m = n == 0 ? 0 : a[0].size();
	assert(ssize(v) == n);
	int cur = 0;
	rep(j, m) {
		int f = -1;
		for (int i = cur; i < n; i++) if (a[i][j] != 0) { f = i; break; }
		if (f == -1) continue;
		if (cur != f) { swap(a[cur], a[f]); swap(v[cur], v[f]); }
		rep(i, n) {
			if (i == cur) continue;
			if (a[i][j] == 0) continue;
			a[i] = a[i] ^ a[cur];
			v[i] ^= v[cur];
		}
		cur++;
		if (cur == n) break;
	}
	BS res(m);
	vector<int> idx(m, -1);
	rep(i, n) {
		int f = -1;
		rep(j, m) if (a[i][j] != 0) { f = j; break; }
		if (f == -1 && v[i] != 0) { return { -1,{} }; }
		if (f != -1) { if (v[i]) res.set(f); idx[f] = i; }
	}
	if (unique) return { 0, {res} };
	vector<BS> vec = { res };
	int dim = 0;
	rep(j, m) {
		if (idx[j] != -1) continue;
		BS w(m); w.set(j);
		rep(i, m) if (idx[i] != -1) if (a[idx[i]][j]) w.set(i);
		dim++; vec.push_back(w);
	}
	return { dim,vec };
}

// a の部分集合で総 xor = x を探す
template<class T> pair<bool, vector<int>> find_linear_comb(vector<T>& a, T x, bool is_good_basis) {
	int sz = ssize(a);
	if (is_good_basis) {
		vector<int> res;
		rep(i, sz) {
			if (chmin(x, x ^ a[i])) res.push_back(i);
		}
		if (x == 0) return { true, res };
		else return { false, {} };
	}
	else {
		vector<pair<T, ull>> basis; vector<int> elem;
		rep(i, sz) {
			T e = a[i];
			T cur = e; ull s = 1ULL << ssize(basis);
			for (auto [b, t] : basis) {
				if (chmin(cur, cur ^ b)) s ^= t;
			}
			if (cur != 0) {
				basis.push_back({ cur, s });
				elem.push_back(i);
			}
		}
		assert(ssize(basis) <= 64);
		{
			T cur = x; ull s = 0;
			for (auto [b, t] : basis) {
				if (chmin(cur, cur ^ b)) s ^= t;
			}
			if (cur != 0) return { false, {} };
			vector<int> res;
			rep(i, ssize(basis)) if (s >> i & 1) res.push_back(elem[i]);
			return { true, res };
		}
	}
}

vector<ll> intersection_F2(vector<ll> a, vector<ll> b) {
	constexpr int D = 30;

	auto f = [&](vector<ll> v) -> array<ll, D> {
		array<ll, D> res{};
		for (ll e : v) {
			for (int i = D - 1; i >= 0; i--) chmin(e, e ^ res[i]);
			if (e) {
				ll t = countr_zero<ull>(bit_floor<ull>(e));
				res[t] = e;
			}
		}
		return res;
	};

	auto fa = f(a), fb = f(b);
	vector<ll> res;
	vector<ll> vec, va;
	rep(i, D) {
		if (fa[i] != 0 && fb[i] != 0) {
			ll target = fa[i] ^ fb[i];
			auto [valid, com] = find_linear_comb(vec, target, true);
			if (valid) {
				ll w = fa[i];
				for (int e : com) w ^= va[e];
				res.push_back(w);
			}
		}
		if (fa[i] != 0) {
			ll w = fa[i];
			int sz = ssize(vec);
			ll from_a = w;
			rep(i, sz) if (chmin(w, w ^ vec[i])) from_a ^= va[i];
			if (w) vec.push_back(w), va.push_back(from_a);
		}
		if (fb[i] != 0) {
			ll w = fb[i];
			int sz = ssize(vec);
			ll from_a = 0;
			rep(i, sz) if (chmin(w, w ^ vec[i])) from_a ^= va[i];
			if (w) vec.push_back(w), va.push_back(from_a);
		}
	}
	return res;
}
