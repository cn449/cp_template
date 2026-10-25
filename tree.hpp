template<class T, auto op, auto e> class Sparsetable {
	int n;
	vector<vector<T>> data;
public:
	Sparsetable() {}
	Sparsetable(vector<T> v) {
		n = ssize(v);
		data.push_back(v);
		for (int i = 2; i <= n; i <<= 1) {
			vector<T> w(n - i + 1);
			for (int j = 0; j <= n - i; j++) {
				w[j] = op(data.back()[j], data.back()[j + i / 2]);
			}
			data.push_back(w);
		}
	}
	T prod(int l, int r) {
		if (l >= r) return e();
		int u = bit_floor<uint>(r - l);
		int d = countr_zero<uint>(u);
		return op(data[d][l], data[d][r - u]);
	}
};

class Tree {
	bool have_doubling_table = false;
	bool have_subtree_size_table = false;
	bool have_heavy_child_table = false;
	bool use_sparse_table_for_lca = false;
	bool have_sparse_table_for_lca = false;
	vector<int> sub_sz, heavy_child;
	static P<int> op_min_pair(P<int> a, P<int> b) { return min(a, b); }
	static P<int> e_min_pair() { return { ten(9), ten(9) }; }
	using sparsetable = Sparsetable<P<int>, op_min_pair, e_min_pair>;
	sparsetable spt;
	void pre_dfs(int v, int p, int d) {
		tour.push_back(v);
		par[v] = p; dep[v] = d;
		for (int i : G[v]) {
			if (i == p) continue;
			pre_dfs(i, v, d + 1);
			tour.push_back(v);
		}
	}
	void init() {
		pre_dfs(root, -1, 0);
		in.resize(n); out.resize(n);
		rep(i, n) in[i] = -1;
		rep(i, 2 * n - 1) {
			if (in[tour[i]] == -1) in[tour[i]] = i;
			out[tour[i]] = i;
		}
	}
public:
	int n, root;
	vector<vector<int>> G;
	vector<int> par, tour, in, out, dep;
	vector<vector<int>> db;
	Tree() {}
	Tree(int n_) {
		n = n_;
		G.resize(n); par.resize(n); dep.resize(n);
	}
	Tree(vector<vector<int>> G_, int root_ = 0) {
		G = move(G_); root = root_;
		n = ssize(G);
		par.resize(n); dep.resize(n);
		init();
	}
	void input_and_init(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		rep(_, n - 1) {
			int u, v;
			cin >> u >> v;
			u--; v--;
			G[u].push_back(v); G[v].push_back(u);
		}
		init();
	}
	void make_random(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		rep(i, n - 1) {
			int u = rng(0, i), v = i + 1;
			G[u].push_back(v); G[v].push_back(u);
		}
		init();
	}
	void make_path(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		int w = root;
		rep(i, n - 1) {
			int u = w, v = w + 1;
			if (v == n) v = 0;
			w++;
			if (w == n) w = 0;
			G[u].push_back(v);
			G[v].push_back(u);
		}
		init();
	}
	void make_almost_path(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		int p = n - n / 100, q = n / 100;
		int w = root;
		rep(i, p - 1) {
			int u = w, v = w + 1;
			if (v == n) v = 0;
			w++;
			if (w == n) w = 0;
			G[u].push_back(v);
			G[v].push_back(u);
		}
		rep(_, q) {
			w++;
			if (w == n) w = 0;
			int u = root + rng(0, p - 1);
			if (u >= n) u -= n;
			G[u].push_back(w);
			G[w].push_back(u);
		}
		init();
	}
	void make_star(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		rep(i, n) {
			int u = root, v = i;
			if (u == v) continue;
			G[u].push_back(v); G[v].push_back(u);
		}
		init();
	}
	void make_almost_star(int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		vector<bool> f(n);
		f[root] = true;
		rep(i, n) {
			int u = root, v = i;
			if (u == v) continue;
			if (rng(0, 99) == 0) continue;
			f[v] = true;
			G[u].push_back(v); G[v].push_back(u);
		}
		rep(i, n) {
			if (!f[i]) {
				int j;
				do { j = rng(0, n - 1); } while (!f[j]);
				G[i].push_back(j);
				G[j].push_back(i);
				f[i] = true;
			}
		}
		init();
	}
	void make_k_ary(int k, int root_ = 0) {
		root = root_;
		assert(0 <= root && root < n);
		for (int i = 1; i < n; i++) {
			int u = (i - 1) / k;
			int v = i;
			u += root; v += root;
			if (u >= n) u -= n;
			if (v >= n) v -= n;
			G[u].push_back(v);
			G[v].push_back(u);
		}
		init();
	}
	// v が u を根とする部分木にありますか？
	bool subtree(int u, int v) {
		assert(0 <= u && u < n); assert(0 <= v && v < n);
		return in[u] <= in[v] && out[v] <= out[u];
	}
	// 長さ 2n, 辺の並びを作る
	vector<P<int>> build_edge() {
		vector<P<int>> edge(2 * n);
		edge[0] = { root, 1 };
		rep(i, 2 * n - 2) {
			int u = tour[i], v = tour[i + 1];
			if (dep[u] < dep[v]) edge[i + 1] = { v, 1 };
			else edge[i + 1] = { u, -1 };
		}
		edge.back() = { root, -1 };
		return edge;
	}
	void build_doubling_table() {
		if (have_doubling_table) return;
		have_doubling_table = true;
		int maxd = *max_element(all(dep));
		int sz = countr_zero<uint>(bit_ceil<uint>(maxd + 2));
		db.resize(sz);
		rep(i, sz) db[i].resize(n + 1);
		rep(j, n) db[0][j] = par[j] == -1 ? n : par[j];
		db[0][n] = n;
		rep(i, sz - 1) rep(j, n + 1) db[i + 1][j] = db[i][db[i][j]];
	}
	void build_sparse_table() {
		if (have_sparse_table_for_lca) return;
		have_sparse_table_for_lca = true;
		vector<P<int>> vec(2 * n - 1);
		rep(i, 2 * n - 1) vec[i] = { dep[tour[i]], tour[i] };
		spt = sparsetable(vec);
	}
	int get_ancestor(int v, int k) {
		assert(0 <= v && v < n); assert(k >= 0);
		if (dep[v] < k) return -1;
		if (!have_doubling_table) build_doubling_table();
		int sz = ssize(db);
		for (int i = sz - 1; i >= 0; i--) if (k >> i & 1) v = db[i][v];
		if (v == n) v = -1;
		return v;
	}
	int lca(int u, int v) {
		assert(0 <= u && u < n); assert(0 <= v && v < n);
		if (use_sparse_table_for_lca) {
			if (!have_sparse_table_for_lca) build_sparse_table();
			int L = in[u], R = in[v];
			if (L > R) swap(L, R);
			return spt.prod(L, R + 1).second;
		}
		else {
			if (!have_doubling_table) build_doubling_table();
			if (dep[u] > dep[v]) swap(u, v);
			int d = dep[v] - dep[u];
			int sz = ssize(db);
			rep(i, sz) if (d >> i & 1) v = db[i][v];
			for (int i = sz - 1; i >= 0; i--) {
				int x = db[i][u], y = db[i][v];
				if (x != y) u = x, v = y;
			}
			if (u != v) u = par[u];
			return u;
		}
	}
	// u から v 方向へ d jump, dist(u, v) < d のとき -1
	int jump(int u, int v, int d) {
		assert(0 <= u && u < n); assert(0 <= v && v < n); assert(d >= 0);
		int l = lca(u, v);
		int d0 = dep[u] - dep[l], d1 = dep[v] - dep[l];
		if (d0 + d1 < d) return -1;
		if (d <= d0) return get_ancestor(u, d);
		else return get_ancestor(v, d0 + d1 - d);
	}
	int dist(int u, int v) {
		assert(0 <= u && u < n); assert(0 <= v && v < n);
		int l = lca(u, v);
		return dep[u] + dep[v] - 2 * dep[l];
	}
	vector<int> diameter() {
		vector<int> D(n), pre(n);
		auto dfs = [&](auto dfs, int v, int p, int d) -> void {
			D[v] = d; pre[v] = p;
			for (int i : G[v]) {
				if (i == p) continue;
				dfs(dfs, i, v, d + 1);
			}
		};
		dfs(dfs, 0, -1, 0);
		int v = max_element(all(D)) - D.begin();
		dfs(dfs, v, -1, 0);
		int u = max_element(all(D)) - D.begin();
		vector<int> res;
		while (u != -1) {
			res.push_back(u);
			u = pre[u];
		}
		return res;
	}
	void calc_subtree_size() {
		if (have_subtree_size_table) return;
		have_subtree_size_table = true;
		sub_sz = vector<int>(n, 1);
		auto dfs = [&](auto dfs, int v, int p) -> void {
			for (int i : G[v]) {
				if (i == p) continue;
				dfs(dfs, i, v);
				sub_sz[v] += sub_sz[i];
			}
		};
		dfs(dfs, root, -1);
	}
	int size(int v) {
		assert(0 <= v && v < n);
		if (!have_subtree_size_table) calc_subtree_size();
		return sub_sz[v];
	}
	void calc_heavy_child() {
		if (have_heavy_child_table) return;
		if (!have_subtree_size_table) calc_subtree_size();
		have_heavy_child_table = true;
		heavy_child = vector<int>(n, -1);
		auto dfs = [&](auto dfs, int v, int p) -> void {
			int mx = -1;
			for (int i : G[v]) {
				if (i == p) continue;
				dfs(dfs, i, v);
				if (chmax(mx, sub_sz[i])) heavy_child[v] = i;
			}
		};
		dfs(dfs, root, -1);
	}
	int heavy(int v) {
		assert(0 <= v && v < n);
		if (!have_heavy_child_table) calc_heavy_child();
		return heavy_child[v];
	}
	vector<P<int>> enum_edge() {
		vector<P<int>> edges;
		rep(i, n) for (int j : G[i]) if (i < j) edges.push_back({ i, j });
		return edges;
	}
	void print() {
		auto edges = enum_edge();
		string s = "";
		for (auto [u, v] : edges) { s += '('; s += to_string(u); s += ','; s += to_string(v); s += ')'; s += ','; }
		if (!s.empty()) s.pop_back();
		cerr << s << '\n';
	}
};