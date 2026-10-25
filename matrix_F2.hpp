class Bitset {
	int n, sz;
	vector<ull> data;
public:
	Bitset(int n_) { n = n_; sz = (n + 63) / 64; data.resize(sz); }
	void set(int pos) { data[pos >> 6] |= 1ULL << (pos & 63); }
	void reset(int pos) { data[pos >> 6] &= ~(1ULL << (pos & 63)); }
	void flip(int pos) { data[pos >> 6] ^= 1ULL << (pos & 63); }
	bool get(int pos) const { return data[pos >> 6] >> (pos & 63) & 1; }
	string to_string() const { string res = ""; rep(i, n) res += get(i) ? '1' : '0'; return res; }
	int count() const { int res = 0; rep(i, sz) res += popcount(data[i]); return res; }
	int size() const { return n; }
	bool operator == (const Bitset& b) const { return data == b.data; }
	bool operator != (const Bitset& b) const { return data != b.data; }
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
	bool operator [](int pos) const { return get(pos); }
};

using BS = Bitset;

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
template<class T> pair<bool, vector<int>> find_linear_comb(vector<T> &a, T x, bool is_good_basis) {
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