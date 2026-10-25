vector<int> bfs(int s, vector<vector<int>>& G) {
	int n = ssize(G);
	queue<int> que;
	que.push(s);
	vector<int> d(n, inf);
	d[s] = 0;
	while (!que.empty()) {
		int q = que.front(); que.pop();
		for (int v : G[q]) {
			if (d[v] == inf) {
				d[v] = d[q] + 1;
				que.push(v);
			}
		}
	}
	return d;
}

vector<ll> dijkstra(int s, vector<vector<P<ll>>>& G) {
	int n = ssize(G);
	vector<ll> d(n, INF);
	d[s] = 0;
	pri_s<P<ll>> pq;
	pq.push({ 0,s });
	while (!pq.empty()) {
		auto [cost, v] = pq.top(); pq.pop();
		if (d[v] < cost) continue;
		for (auto [nv, ncost] : G[v]) if (chmin(d[nv], d[v] + ncost)) pq.push({ d[nv], nv });
	}
	return d;
}

// 負閉路なし
vector<vector<ll>> APSP(vector<vector<P<ll>>> G) {
	int n = ssize(G);
	vector<ll> pot(n);
	rep(_, n) {
		bool upd = false;
		rep(v, n) for (auto& [u, w] : G[v]) {
			if (chmin(pot[u], pot[v] + w)) upd = true;
		}
		if (!upd) break;
	}
	rep(v, n) for (auto& [u, w] : G[v]) w += pot[v] - pot[u];
	vector d(n, vector<ll>());
	rep(i, n) d[i] = dijkstra(i, G);
	rep(i, n) rep(j, n) d[i][j] += pot[j] - pot[i];
	return d;
}

// 存在するか, （存在するならば）サイクルを 1 つ取る
// CF905-E, ARC079-F
pair<bool, vector<int>> find_cycle(vector<vector<int>> G, bool directed) {
	int n = ssize(G);
	rep(i, n) {
		sort(all(G[i]));
		int sz = ssize(G[i]);
		// 自己ループ
		for (int v : G[i]) if (i == v) return { true, { int(i) } };
		// 多重辺
		if (!directed) rep(j, sz - 1) if (G[i][j] == G[i][j + 1]) return { true, { int(i), G[i][j] } };
	}
	vector<int> state(n);
	vector<int> res;
	vector<int> pre(n, -1);
	int start = -1, goal = -1;

	auto dfs = [&](auto dfs, int v, int p) -> void {
		state[v] = 1;
		pre[v] = p;
		for (int nxt : G[v]) {
			if (start != -1) return;
			if (!directed && nxt == p) continue;
			if (state[nxt] == 0) dfs(dfs, nxt, v);
			else if (state[nxt] == 1) {
				start = v;
				goal = nxt;
				return;
			}
		}
		state[v] = 2;
	};

	rep(i, n) if (state[i] == 0 && start == -1) dfs(dfs, i, -1);
	if (start == -1) return { false, {} };
	int cur = start;
	while (true) {
		res.push_back(cur);
		if (cur == goal) break;
		cur = pre[cur];
	}
	if (directed) reverse(all(res));
	return { true, res };
}

// 存在するか, 頂点列, 辺の列
tuple<bool, vector<int>, vector<int>> eulerian_path(vector<P<int>> E, int s, bool directed) {
	int n = 0;
	for (auto& e : E) chmax(n, max(e.first, e.second));
	n++;
	vector G(n, vector<P<int>>());
	for (auto [u, v] : E) {
		int sz_v = ssize(G[v]);
		G[u].push_back({ v, sz_v });
		if (!directed) {
			int sz_u = ssize(G[u]) - 1;
			G[v].push_back({ u, sz_u });
		}
	}
	vector<int> v_ord, tmp;
	tmp.push_back(s);
	while (!tmp.empty()) {
		int u = tmp.back();
		if (G[u].empty()) {
			v_ord.push_back(u);
			tmp.pop_back();
		}
		else {
			auto [v, idx] = G[u].back();
			G[u].pop_back();
			if (idx == -1) continue;
			if (!directed) G[v][idx].second = -1;
			tmp.push_back(v);
		}
	}
	reverse(all(v_ord));
	int sz = ssize(v_ord);
	if (sz != ssize(E) + 1) return { false, {},{} };
	map<P<int>, vector<int>> mp;
	rep(i, ssize(E)) {
		P<int> p = E[i];
		if (!directed && p.first > p.second) swap(p.first, p.second);
		mp[p].push_back(i);
	}
	vector<int> e_ord;
	rep(i, sz - 1) {
		int u = v_ord[i], v = v_ord[i + 1];
		if (!directed && u > v) swap(u, v);
		if (mp[{u, v}].empty()) return { false, {},{} };
		int t = mp[{u, v}].back(); mp[{u, v}].pop_back();
		e_ord.push_back(t);
	}
	return { true, v_ord, e_ord };
}

vector<int> topological_sort(vector<vector<int>>& G) {
	int n = ssize(G);
	vector<int> indegree(n);
	rep(i, n) for (int j : G[i]) indegree[j]++;
	vector<int> res;
	queue<int> que;
	for (int i = 0; i < n; i++) if (indegree[i] == 0) que.push(i);
	while (!que.empty()) {
		int ver = que.front(); que.pop();
		res.push_back(ver);
		for (int i : G[ver]) {
			indegree[i]--;
			if (indegree[i] == 0) que.push(i);
		}
	}
	return res;
}