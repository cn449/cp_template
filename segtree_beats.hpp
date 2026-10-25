// ACL based
// verify : https://judge.yosupo.jp/problem/range_chmin_chmax_add_range_sum
// submission : https://judge.yosupo.jp/submission/407066
template <class S, auto op, auto e, class F, auto mapping, auto composition, auto id> struct segtree_beats {
	static_assert(is_convertible_v<decltype(op), function<S(S, S)>>, "op must work as S(S, S)");
	static_assert(is_convertible_v<decltype(e), function<S()>>, "e must work as S()");
	static_assert(is_convertible_v<decltype(mapping), function<S(F, S)>>, "mapping must work as F(F, S)");
	static_assert(is_convertible_v<decltype(composition), function<F(F, F)>>, "compostiion must work as F(F, F)");
	static_assert(is_convertible_v<decltype(id), function<F()>>, "id must work as F()");
public:
	segtree_beats() : segtree_beats(0) {}
	explicit segtree_beats(int n) : segtree_beats(vector<S>(n, e())) {}
	explicit segtree_beats(const vector<S>& v) : _n(ssize(v)) {
		size = bit_ceil<uint>(_n);
		log = countr_zero<uint>(size);
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
		if (d[k].fail) { push(k); update(k); }
	}
	void push(int k) {
		all_apply(2 * k, lz[k]);
		all_apply(2 * k + 1, lz[k]);
		lz[k] = id();
	}
};

struct S {
	ll sum, mn, mn2, mncnt, mx, mx2, mxcnt, len;
	bool fail;
};

S op(S a, S b) {
	S res;
	res.sum = a.sum + b.sum;
	res.mn = min(a.mn, b.mn);
	res.mx = max(a.mx, b.mx);
	if (a.mn < b.mn) {
		res.mn2 = min(a.mn2, b.mn);
		res.mncnt = a.mncnt;
	}
	else if (a.mn == b.mn) {
		res.mn2 = min(a.mn2, b.mn2);
		res.mncnt = a.mncnt + b.mncnt;
	}
	else {
		res.mn2 = min(a.mn, b.mn2);
		res.mncnt = b.mncnt;
	}
	if (a.mx > b.mx) {
		res.mx2 = max(a.mx2, b.mx);
		res.mxcnt = a.mxcnt;
	}
	else if (a.mx == b.mx) {
		res.mx2 = max(a.mx2, b.mx2);
		res.mxcnt = a.mxcnt + b.mxcnt;
	}
	else {
		res.mx2 = max(a.mx, b.mx2);
		res.mxcnt = b.mxcnt;
	}
	res.len = a.len + b.len;
	res.fail = false;
	return res;
}

S e() {
	return { 0, INF, INF, 0, -INF, -INF, 0, 0, false };
}

struct F { ll add, lb, ub; }; // x -> clamp(x + add, lb, ub), lb <= ub

S mapp(F f, S a) {
	a.sum += a.len * f.add;
	a.mn += f.add;
	a.mn2 += f.add;
	a.mx += f.add;
	a.mx2 += f.add;
	if (a.mn < f.lb && f.lb < a.mn2) {
		a.sum += (f.lb - a.mn) * a.mncnt;
		if (a.mx == a.mn) a.mx = f.lb;
		else if (a.mx2 == a.mn) a.mx2 = f.lb;
		a.mn = f.lb;
	}
	else if (a.mn2 <= f.lb) a.fail = true;
	if (a.mx2 < f.ub && f.ub < a.mx) {
		a.sum -= (a.mx - f.ub) * a.mxcnt;
		if (a.mn == a.mx) a.mn = f.ub;
		else if (a.mn2 == a.mx) a.mn2 = f.ub;
		a.mx = f.ub;
	}
	else if (f.ub <= a.mx2) a.fail = true;
	return a;
}

F comp(F g, F f) {
	f.add += g.add;
	f.lb += g.add;
	f.ub += g.add;
	if (g.ub <= f.lb) {
		f.lb = g.ub;
		f.ub = g.ub;
	}
	else if (f.ub <= g.lb) {
		f.lb = g.lb;
		f.ub = g.lb;
	}
	else {
		chmax(f.lb, g.lb);
		chmin(f.ub, g.ub);
	}
	return f;
}

F id() { return { 0, -INF, INF }; }