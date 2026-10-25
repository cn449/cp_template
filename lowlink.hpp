// 多重辺, 非連結に対応
class lowlink {
public:
	int n; // 頂点数
	int component; // 連結成分数
	vector<vector<pair<int, bool>>> G; // 接続する頂点, 橋か
	vector<vector<int>> dfs_tree; // 親 -> 子方向の辺だけ張ってます
	vector<bool> root;
	vector<bool> used;
	vector<int> ord, low, par;
	vector<P<int>> bridge;
	void dfs(int u, int& k, int p) {
		ord[u] = low[u] = k++; par[u] = p;
		used[u] = true;
		rep(i, ssize(G[u])) {
			int v = G[u][i].first;
			if (!used[v]) {
				dfs(v, k, u);
				dfs_tree[u].push_back(v);
				chmin(low[u], low[v]);
				if (ord[u] < low[v] && !(i + 1 < ssize(G[u]) && G[u][i].first == G[u][i + 1].first)) bridge.push_back(minmax(u, v));
			}
			else if (v != p) chmin(low[u], ord[v]);
		}
	}
	bool is_bridge(int u, int v) {
		P<int> p = minmax(u, v);
		return binary_search(all(bridge), p);
	}
	lowlink(vector<vector<int>> g) {
		n = ssize(g);
		G.resize(n);
		rep(i, n) {
			for (int j : g[i]) G[i].push_back({ j, false });
			sort(all(G[i]));
		}
		dfs_tree.resize(n);
		root.resize(n); used.resize(n);
		ord.resize(n); low.resize(n); par.resize(n);
		int k = 0;
		rep(i, n) if (!used[i]) {
			root[i] = true;
			component++;
			dfs(i, k, -1);
		}
		sort(all(bridge));
		rep(i, n) for (auto& [j, flag] : G[i]) if (is_bridge(i, j)) flag = true;
	}
	// 頂点 u を削除したとき、連結成分数はいくつになりますか？
	int count_component(int u) {
		int c = 0;
		for (int v : dfs_tree[u]) {
			if (ord[u] <= low[v]) c++;
		}
		return component + c - (root[u] ? 1 : 0);
	}
	bool is_articulation(int u) {
		int c = 0;
		for (int v : dfs_tree[u]) {
			if (ord[u] <= low[v]) c++;
		}
		int w = c - (root[u] ? 1 : 0);
		return w > 0;
	}
	vector<vector<int>> two_edge_cc() {
		vector<bool> f(n);
		vector<vector<int>> res;
		rep(i, n) {
			if (f[i]) continue;
			f[i] = true;
			vector<int> vec;
			vec.push_back(i);
			queue<int> que;
			que.push(i);
			while (!que.empty()) {
				int q = que.front(); que.pop();
				for (auto [v, bad] : G[q]) {
					if (bad or f[v]) continue;
					f[v] = true;
					vec.push_back(v);
					que.push(v);
				}
			}
			res.push_back(vec);
		}
		return res;
	}
	vector<vector<int>> two_vertex_cc() {
		vector<bool> f(n);
		vector<vector<int>> res;
		vector<int> seq(n);
		rep(i, n) seq[i] = i;
		sort(all(seq), [&](int u, int v) { return ord[u] < ord[v]; });
		for (int i : seq) {
			if (par[i] == -1) {
				if (ssize(G[i]) == 0) res.push_back({ i });
				continue;
			}
			if (f[i]) continue;
			vector<int> vec;
			auto rec = [&](auto rec, int v, int p) -> void {
				f[v] = true;
				vec.push_back(v);
				for (int u : dfs_tree[v]) {
					if (low[u] < ord[v]) rec(rec, u, v);
				}
			};
			vec.push_back(par[i]);
			rec(rec, i, par[i]);
			res.push_back(vec);
		}
		return res;
	}
};