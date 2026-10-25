class SCC {
	int n, scc_sz;
	vector<vector<int>> G, RG;
	vector<int> vs;
	vector<bool> used;
	vector<int> cmp;
	bool build = false;
public:
	SCC(int n_) {
		n = n_;
		G = RG = vector<vector<int>>(n, vector<int>());
		cmp = vector<int>(n);
		used = vector<bool>(n);
	}
	void add_edge(int from, int to) {
		G[from].push_back(to);
		RG[to].push_back(from);
	}
	void dfs(int v) {
		used[v] = true;
		for (int i : G[v]) if (!used[i]) dfs(i);
		vs.push_back(v);
	}
	void rdfs(int v, int k) {
		used[v] = true;
		cmp[v] = k;
		for (int i : RG[v]) if (!used[i]) rdfs(i, k);
	}
	void solve() {
		assert(!build);
		build = true;
		rep(v, n) if (!used[v]) dfs(v);
		rep(i, n) used[i] = false;
		int k = 0;
		for (int i = n - 1; i >= 0; i--) if (!used[vs[i]]) rdfs(vs[i], k++);
		scc_sz = k;
	}
	int size() {
		if (!build) solve();
		return scc_sz;
	}
	bool same(int u, int v) {
		if (!build) solve();
		return cmp[u] == cmp[v];
	}
	int get(int v) {
		if (!build) solve();
		return cmp[v];
	}
	// DAG, index small -> large 辺のみになる
	vector<vector<int>> get_graph() {
		if (!build) solve();
		vector H(scc_sz, vector<int>());
		rep(i, n) {
			for (int j : G[i]) {
				int u = cmp[i], v = cmp[j];
				if (u < v) H[u].push_back(v);
			}
		}
		rep(i, scc_sz) {
			sort(all(H[i]));
			H[i].erase(unique(all(H[i])), H[i].end());
		}
		return H;
	}
};