template<class T> class Dinic {
	struct edge {
		int to;
		T cap, icap;
		int rev;
	};
	vector<vector<edge>> graph;
	vector<int> level, iter;
public:
	Dinic(int n) {
		graph = vector<vector<edge>>(n, vector<edge>());
		level = iter = vector<int>(n);
	}
	void add_edge(int from, int to, T cap) {
		graph[from].push_back({ to,cap,cap,(T)graph[to].size() });
		graph[to].push_back({ from,0,0,(T)graph[from].size() - 1 });
	}
	void bfs(int s) {
		for (int& i : level) i = -1;
		queue<int> que;
		level[s] = 0;
		que.push(s);
		while (!que.empty()) {
			int v = que.front(); que.pop();
			for (int i = 0; i < graph[v].size(); i++) {
				edge& e = graph[v][i];
				if (e.cap > 0 && level[e.to] < 0) {
					level[e.to] = level[v] + 1;
					que.push(e.to);
				}
			}
		}
	}
	T dfs(int v, int t, T f) {
		if (v == t) return f;
		for (int& i = iter[v]; i < graph[v].size(); i++) {
			edge& e = graph[v][i];
			if (e.cap > 0 && level[v] < level[e.to]) {
				T d = dfs(e.to, t, min(f, e.cap));
				if (d > 0) {
					e.cap -= d;
					graph[e.to][e.rev].cap += d;
					return d;
				}
			}
		}
		return 0;
	}
	T max_flow(int s, int t) {
		T flow = 0;
		while (true) {
			bfs(s);
			if (level[t] < 0) return flow;
			for (int& i : iter) i = 0;
			T f = 0;
			while ((f = dfs(s, t, inf)) > 0) flow += f;
		}
	}
	vector<P<int>> bipartite_matching(int s, int t) {
		vector<P<int>> res;
		max_flow(s, t);
		for (int i = 0; i < graph.size(); i++) {
			if (i == s) continue;
			for (edge e : graph[i]) {
				if (e.cap == 0 && e.icap == 1 && e.to != t) {
					res.push_back({ i,e.to });
				}
			}
		}
		return res;
	}
};