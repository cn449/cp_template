// static, offline
// T : 重みの型
// 座標の与え方は lrdu 順 (x \in [l, r), y \in [d, u))
template<class T> vector<T> rectangle_sum(vector<tuple<ll, ll, T>> pt, vector<array<ll, 4>> qu) {
	int n = ssize(pt);
	int q = ssize(qu);
	vector<ll> X, Y;
	for (auto [x, y, w] : pt) {
		X.push_back(x); Y.push_back(y);
	}
	for (auto [l, r, d, u] : qu) {
		X.push_back(l); X.push_back(r);
		Y.push_back(d); Y.push_back(u);
	}
	sort(all(X)); X.erase(unique(all(X)), X.end());
	sort(all(Y)); Y.erase(unique(all(Y)), Y.end());
	int sx = ssize(X), sy = ssize(Y);
	vector pts(sx, vector<pair<int, T>>());
	for (auto [x, y, w] : pt) {
		x = lower_bound(all(X), x) - X.begin();
		y = lower_bound(all(Y), y) - Y.begin();
		pts[x].push_back({ y, w });
	}
	vector add(sx, vector<array<int, 4>>());
	rep(i, q) {
		auto [l, r, d, u] = qu[i];
		l = lower_bound(all(X), l) - X.begin();
		r = lower_bound(all(X), r) - X.begin();
		d = lower_bound(all(Y), d) - Y.begin();
		u = lower_bound(all(Y), u) - Y.begin();
		add[l].push_back({ int(i), int(d), int(u), -1 });
		add[r].push_back({ int(i), int(d), int(u), 1 });
	}
	vector<T> ans(q);
	fenwick_tree<ll> F(sy);
	rep(i, sx) {
		for (auto [ind, l, r, coeff] : add[i]) ans[ind] += F.sum(l, r) * coeff;
		for (auto [y, w] : pts[i]) F.add(y, w);
	}
	return ans;
}

struct S { ll mn, mn_len = 0; };
S op(S a, S b) {
	S res;
	res.mn = min(a.mn, b.mn);
	if (res.mn == a.mn) res.mn_len += a.mn_len;
	if (res.mn == b.mn) res.mn_len += b.mn_len;
	return res;
}
S e() { return { inf, 0 }; }
S mapp(ll f, S s) { return { s.mn + f, s.mn_len }; }
ll comp(ll g, ll f) { return g + f; }
ll id() { return 0; }

ll area_of_union_of_rectangles(vector<P<P<ll>>> rects) {
	if (rects.empty()) return 0;
	vector<ll> X, Y;
	for (auto& [x, y] : rects) {
		auto& [xl, xr] = x; X.push_back(xl); X.push_back(xr);
		auto& [yl, yr] = y; Y.push_back(yl); Y.push_back(yr);
	}
	sort(all(X)); X.erase(unique(all(X)), X.end());
	sort(all(Y)); Y.erase(unique(all(Y)), Y.end());
	int n = ssize(X), m = ssize(Y);
	vector<S> init(m - 1);
	rep(i, m - 1) init[i] = { 0, Y[i + 1] - Y[i] };
	lazy_segtree<S, op, e, ll, mapp, comp, id> seg(init);
	vector add(n, vector<array<ll, 3>>());
	for (auto& [x, y] : rects) {
		auto& [xl, xr] = x;
		xl = lower_bound(all(X), xl) - X.begin();
		xr = lower_bound(all(X), xr) - X.begin();
		auto& [yl, yr] = y;
		yl = lower_bound(all(Y), yl) - Y.begin();
		yr = lower_bound(all(Y), yr) - Y.begin();
		add[xl].push_back({ yl, yr, 1 });
		add[xr].push_back({ yl, yr, -1 });
	}
	ll res = 0;
	ll w = Y.back() - Y[0];
	rep(i, n - 1) {
		for (auto [l, r, c] : add[i]) seg.apply(l, r, c);
		auto [mn, mn_len] = seg.all_prod();
		ll non_zero = mn == 0 ? (w - mn_len) : w;
		res += non_zero * (X[i + 1] - X[i]);
	}
	return res;
}