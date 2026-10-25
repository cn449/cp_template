class mincostflow {
	struct edge {
		ll to, cap, cost, rev;
	};
	vector<vector<edge>> G;
	vector<ll> h, dist, prevv, preve;
	bool negative_edge = false;
public:
	mincostflow(int n) {
		G = vector<vector<edge>>(n, vector<edge>());
		h = dist = prevv = preve = vector<ll>(n);
	}
	void add_edge(int from, int to, ll cap, ll cost) {
		G[from].push_back({ to, cap, cost, ssize(G[to]) });
		G[to].push_back({ from, 0, -cost, ssize(G[from]) - 1 });
		if (cost < 0) negative_edge = true;
	}
	P<ll> flow_dijkstra(int s, int t, ll f) {
		// flow, cost を返す
		// 負辺なし、O(FElogV)
		ll flow = f, cost = 0;
		while (f > 0) {
			for (ll& e : dist) e = INF;
			dist[s] = 0;
			pri_s<P<ll>> que;
			que.push({ 0,s });
			while (!que.empty()) {
				auto p = que.top(); que.pop();
				int v = p.second;
				if (dist[v] < p.first) continue;
				for (int i = 0; i < ssize(G[v]); i++) {
					edge& e = G[v][i];
					if (e.cap > 0 && dist[e.to] > dist[v] + e.cost + h[v] - h[e.to]) {
						dist[e.to] = dist[v] + e.cost + h[v] - h[e.to];
						prevv[e.to] = v, preve[e.to] = i;
						que.push({ dist[e.to],e.to });
					}
				}
			}
			if (dist[t] == INF) return { flow - f, cost };
			rep(v, ssize(h)) h[v] += dist[v];
			ll d = f;
			for (int v = t; v != s; v = prevv[v]) chmin(d, G[prevv[v]][preve[v]].cap);
			f -= d;
			cost += d * h[t];
			for (int v = t; v != s; v = prevv[v]) {
				edge& e = G[prevv[v]][preve[v]];
				e.cap -= d;
				G[v][e.rev].cap += d;
			}
		}
		return { flow, cost };
	}
	P<ll> flow_bellman_ford(int s, int t, ll f) {
		// flow, cost を返す
		// 負辺あり、負閉路なし O(FVE)
		ll flow = f, cost = 0;
		while (f > 0) {
			for (ll& e : dist) e = INF;
			dist[s] = 0;
			bool update = true;
			while (update) {
				update = false;
				for (int v = 0; v < ssize(G); v++) {
					if (dist[v] == INF) continue;
					for (int i = 0; i < ssize(G[v]); i++) {
						edge& e = G[v][i];
						if (e.cap > 0 && dist[e.to] > dist[v] + e.cost) {
							dist[e.to] = dist[v] + e.cost;
							prevv[e.to] = v, preve[e.to] = i;
							update = true;
						}
					}
				}
			}
			if (dist[t] == INF) return { flow - f, cost };
			ll d = f;
			for (int v = t; v != s; v = prevv[v]) chmin(d, G[prevv[v]][preve[v]].cap);
			f -= d;
			cost += d * dist[t];
			for (int v = t; v != s; v = prevv[v]) {
				edge& e = G[prevv[v]][preve[v]];
				e.cap -= d;
				G[v][e.rev].cap += d;
			}
		}
		return { flow, cost };
	}
	P<ll> flow_dag(int s, int t, ll f = INF) {
		// flow, cost を返す
		// DAG O(FElogV)
		ll flow = f, cost = 0;
		bool first_step = true;
		while (f > 0) {
			for (ll& e : dist) e = INF;
			if (first_step) {
				int n = ssize(G);
				vector<int> deg(n);
				rep(i, n) for (auto [to, cap, cost, _] : G[i]) if (cap > 0) deg[to]++;
				dist[s] = 0;
				queue<int> que;
				rep(i, n) if (deg[i] == 0) que.push(i);
				while (!que.empty()) {
					int v = que.front(); que.pop();
					for (int i = 0; i < ssize(G[v]); i++) {
						edge& e = G[v][i];
						if (e.cap == 0) continue;
						if (dist[e.to] > dist[v] + e.cost + h[v] - h[e.to]) {
							dist[e.to] = dist[v] + e.cost + h[v] - h[e.to];
							prevv[e.to] = v, preve[e.to] = i;
						}
						deg[e.to]--;
						if (deg[e.to] == 0) que.push(e.to);
					}
				}
				first_step = false;
			}
			else {
				dist[s] = 0;
				pri_s<P<ll>> que;
				que.push({ 0,s });
				while (!que.empty()) {
					auto p = que.top(); que.pop();
					int v = p.second;
					if (dist[v] < p.first) continue;
					for (int i = 0; i < ssize(G[v]); i++) {
						edge& e = G[v][i];
						if (e.cap > 0 && dist[e.to] > dist[v] + e.cost + h[v] - h[e.to]) {
							dist[e.to] = dist[v] + e.cost + h[v] - h[e.to];
							prevv[e.to] = v, preve[e.to] = i;
							que.push({ dist[e.to],e.to });
						}
					}
				}
			}
			if (dist[t] == INF) return { flow - f, cost };
			rep(v, ssize(h)) h[v] += dist[v];
			ll d = f;
			for (int v = t; v != s; v = prevv[v]) chmin(d, G[prevv[v]][preve[v]].cap);
			f -= d;
			cost += d * h[t];
			for (int v = t; v != s; v = prevv[v]) {
				edge& e = G[prevv[v]][preve[v]];
				e.cap -= d;
				G[v][e.rev].cap += d;
			}
		}
		return { flow, cost };
	}
	P<ll> flow(int s, int t, ll f = INF) { return negative_edge ? flow_bellman_ford(s, t, f) : flow_dijkstra(s, t, f); }
	vector<array<ll, 3>> edges() {
		vector<array<ll, 3>> res;
		for (int i = 0; i < ssize(G); i++) {
			for (edge e : G[i]) {
				if (e.cap > 0) res.push_back({ i,e.to,e.cap });
			}
		}
		return res;
	}
};