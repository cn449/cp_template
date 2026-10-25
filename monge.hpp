// min を達成する index（広義単調増加）を返す
// monotone であれば totally monotone でなくても ok
template<class F> vector<int> monotone_minima(int n, int m, F f) {
	using val_type = invoke_result_t<F&, int, int>;
	vector<int> res(n);

	auto rec = [&](auto rec, int L, int R, int pl, int pr) -> void {
		if (L == R) return;
		int M = midpoint(L, R);
		val_type mn = numeric_limits<val_type>::max();
		For(j, pl, pr) if (chmin(mn, f(M, j))) res[M] = j;
		rec(rec, L, M, pl, res[M] + 1);
		rec(rec, M + 1, R, res[M], pr);
	};

	rec(rec, 0, n, 0, m);
	return res;
}

template<class T> vector<T> convex_arbitrary_minplus(vector<T>& a, vector<T>& b) {
	int n = ssize(a), m = ssize(b);

	auto f = [&](int i, int j) -> T {
		int x = i - j, y = j;
		if (0 <= x && x < n && 0 <= y && y < m) return a[x] + b[y];
		else return numeric_limits<T>::max();
	};

	auto ind = monotone_minima(n + m - 1, m, f);
	vector<T> ans(n + m - 1);
	rep(i, n + m - 1) ans[i] = f(i, ind[i]);
	return ans;
}

// concave をたくさん畳み込んで何かを求めるというときには CF681-D のような方法が使えることもある
// 調和級数だったり分割統治だったり

template<class T> vector<T> concave_arbitrary_minplus_inner(vector<T>& a, vector<T>& b) {
	int n = ssize(a), m = ssize(b);
	assert(n >= m);
	vector<T> ans(n + m - 1, numeric_limits<T>::max());

	auto f = [&](int i, int j) -> T {
		j = m - 1 - j;
		int x = i - j, y = j;
		if (0 <= x && x < n && 0 <= y && y < m) return a[x] + b[y];
		else return numeric_limits<T>::max();
	};

	auto rect = [&](int lx, int rx, int ly, int ry) {
		auto g = [&](int i, int j) -> T { return f(lx + i, ly + j); };
		auto ind = monotone_minima(rx - lx, ry - ly, g);
		rep(i, rx - lx) chmin(ans[i + lx], g(i, ind[i]));
	};

	auto lower = [&](auto lower, int lx, int rx, int ly, int ry) {
		assert(rx - lx == ry - ly);
		if (lx + 1 == rx) {
			chmin(ans[lx], f(lx, ly));
			return;
		}
		int mx = midpoint(lx, rx), my = ry - mx + lx;
		rect(mx, rx, my, ry);
		lower(lower, lx, mx, my, ry);
		lower(lower, mx, rx, ly, my);
	};

	auto upper = [&](auto upper, int lx, int rx, int ly, int ry) {
		assert(rx - lx == ry - ly);
		if (lx + 1 == rx) {
			chmin(ans[lx], f(lx, ly));
			return;
		}
		int mx = midpoint(lx, rx), my = ry - mx + lx;
		rect(lx, mx, ly, my);
		upper(upper, lx, mx, my, ry);
		upper(upper, mx, rx, ly, my);
	};

	lower(lower, 0, m, 0, m);
	if (m <= n - 1) rect(m, n - 1, 0, m);
	upper(upper, n - 1, n + m - 1, 0, m);

	return ans;
}

template<class T> vector<T> concave_arbitrary_minplus(vector<T>& a, vector<T>& b) {
	int n = ssize(a), m = ssize(b);
	vector<T> ans(n + m - 1, numeric_limits<T>::max());
	for (int L = 0; L < m; L += n) {
		int R = min(L + n, m);
		vector<T> c(R - L);
		rep(i, R - L) c[i] = b[i + L];
		auto res = concave_arbitrary_minplus_inner(a, c);
		rep(i, ssize(res)) chmin(ans[L + i], res[i]);
	}
	return ans;
}

template<class F> vector<invoke_result_t<F&, int, int>> monge_dp_slow(int n, F f) {
	using val_type = invoke_result_t<F&, int, int>;
	vector<val_type> dp(n + 1);
	For(i, 1, n + 1) dp[i] = numeric_limits<val_type>::max();

	auto rec = [&](auto rec, int L, int R) -> void {
		if (L + 1 == R) return;
		int M = midpoint(L, R);
		rec(rec, L, M);
		auto g = [&](int i, int j) -> val_type { return dp[L + j] + f(L + j, M + i); };
		auto ind = monotone_minima(R - M, M - L, g);
		rep(i, R - M) chmin(dp[M + i], g(i, ind[i]));
		rec(rec, M, R);
	};

	rec(rec, 0, n + 1);
	return dp;
}

template<class F> vector<invoke_result_t<F&, int, int>> monge_dp(int n, F f) {
	using val_type = invoke_result_t<F&, int, int>;
	vector<val_type> dp(n + 1);
	For(i, 1, n + 1) dp[i] = numeric_limits<val_type>::max();
	vector<int> ind(n + 1);

	auto upd = [&](int i, int j) {
		if (i <= j) return;
		val_type nxt = dp[j] + f(j, i);
		if (chmin(dp[i], nxt)) ind[i] = j;
	};

	auto rec = [&](auto rec, int L, int R) -> void {
		if (L + 1 == R) return;
		int M = midpoint(L, R);
		For(i, ind[L], ind[R] + 1) upd(M, i);
		rec(rec, L, M);
		For(i, L + 1, M + 1) upd(R, i);
		rec(rec, M, R);
	};

	upd(n, 0);
	rec(rec, 0, n);
	return dp;
}

// 1 個辺を辿るたびにコストが t 加算されると思って解き、td を引く
// その値を f(t) とすると、定め方から明らかに max f(t) <= 求める解
// g(d) : d 辺使ったときのコストとして、g が下に凸なら max f(t) = 求める解になる
// コスト関数が monge のとき g は下に凸になる
// g(d-1) + g(d+1) >= g(d) + g(d) を示す
// d-1, d+1 辺辿るときの最適解を取って、そのときの経由点を A_0,...,A_{d-1}, B_0,...,B_{d+1} とする
// A_i, B_{i+1}, B_{i+2}, A_{i+1} がこの順で現れるところがあるので、A_i -> B_{i+2}, B_{i+1} -> A_{i+1} と組み替えれば monge 性より ok
// g(d) - g(d+1) <= x <= g(d-1) - g(d) の範囲にある x が最適なペナルティなので、この範囲を絶対含むように上下限を設定して探索する
// 辺コストの絶対値の最大値を C として、1 本増やしたり減らしたりしても 3C までしか改悪にならないから [-3C, 3C] とか
// 特にすべてが非負のときは辺を取り除いてコストが増えないので [-2C, C] でよい（が、上の評価から計算量的な改善には全然なってない）

// ABC355-G, CF901-D
template<class F> invoke_result_t<F&, int, int> monge_d_edge_dp(int n, F f, invoke_result_t<F&, int, int> lb, invoke_result_t<F&, int, int> ub, int d) {
	using val_type = invoke_result_t<F&, int, int>;
	auto get = [&](val_type t) { return monge_dp(n, [&](int i, int j) {return f(i, j) + t; }).back() - t * d;  };
	if (is_integral_v<val_type>) {
		auto safe_div3 = [&](val_type v) {
			if (v >= 0) return v / 3;
			else return (v - 2) / 3;
		};
		while (ub - lb > 2) {
			val_type m1 = safe_div3(lb + lb + ub);
			val_type m2 = safe_div3(lb + ub + ub);
			val_type f1 = get(m1), f2 = get(m2);
			if (f1 >= f2) ub = m2;
			else lb = m1;
		}
		return max({ get(lb), get(lb + 1), get(ub) });
	}
	else {
		rep(_, 100) {
			val_type m1 = (lb + lb + ub) / 3;
			val_type m2 = (lb + ub + ub) / 3;
			val_type f1 = get(m1), f2 = get(m2);
			if (f1 >= f2) ub = m2;
			else lb = m1;
		}
		return get(lb);
	}
}

// これ手元だと 2 次元 vector の宣言を単に vector にするとコンパイル通らなくて謎
// 要検証だとは思うけど少し調べただけでは原因不明
template<class F> ll knuth_yao_speedup(int n, F f) {
	vector<vector<ll>> dp(n + 1, vector<ll>(n + 1, INF));
	vector<vector<int>> opt(n + 1, vector<int>(n + 1));
	rep(i, n) {
		dp[i][i + 1] = 0;
		opt[i][i + 1] = i;
	}
	for (int l = n - 2; l >= 0; l--) {
		for (int r = l + 2; r <= n; r++) {
			int L = opt[l][r - 1];
			int R = min(r, opt[l + 1][r] + 1);
			For(m, L, R) {
				ll w = dp[l][m] + dp[m][r] + f(l, r);
				if (chmin(dp[l][r], w)) opt[l][r] = m;
			}
		}
	}
	return dp[0][n];
}