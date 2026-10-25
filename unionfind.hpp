class unionfind {
	vector<int> par;
	vector<int> sz;
public:
	unionfind(int n) {
		par = vector<int>(n);
		for (int i = 0; i < n; i++) par[i] = i;
		sz = vector<int>(n, 1);
	}
	unionfind(int n, vector<int> a) {
		par = vector<int>(n);
		for (int i = 0; i < n; i++) par[i] = i;
		sz = a;
	}
	int find(int x) {
		if (par[x] == x) return x;
		else return par[x] = find(par[x]);
	}
	int size(int x) { return sz[find(x)]; }
	bool same(int x, int y) { return find(x) == find(y); }
	void unite(int x, int y) {
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (sz[x] < sz[y]) {
			par[x] = y;
			sz[y] += sz[x];
		}
		else {
			par[y] = x;
			sz[x] += sz[y];
		}
	}
};

// CF1094F, ECR178-G

template<class T> class undoable_vector {
	vector<T> data;
	vector<pair<int, T>> hist;
public:
	undoable_vector() {}
	undoable_vector(int n) { data.resize(n); }
	void set(int pos, T val) {
		hist.push_back({ pos, data[pos] });
		data[pos] = val;
	}
	T get(int pos) { return data[pos]; }
	void undo() {
		auto [pos, val] = hist.back();
		hist.pop_back();
		data[pos] = val;
	}
	T operator [](int i) { return get(i); }
};
 
class undoable_unionfind {
	undoable_vector<int> par;
	undoable_vector<int> sz;
	vector<bool> hist;
	int components;
public:
	undoable_unionfind(int n) {
		par = undoable_vector<int>(n);
		rep(i, n) par.set(i, i);
		sz = undoable_vector<int>(n);
		rep(i, n) sz.set(i, 1);
		components = n;
	}
	int find(int x) {
		if (par[x] == x) return x;
		else return find(par[x]);
	}
	int size(int x) { return sz[find(x)]; }
	bool same(int x, int y) { return find(x) == find(y); }
	void unite(int x, int y) {
		x = find(x); y = find(y);
		if (x == y) {
			hist.push_back(false);
			return;
		}
		components--;
		hist.push_back(true);
		if (sz[x] < sz[y]) {
			par.set(x, y);
			sz.set(y, sz[x] + sz[y]);
		}
		else {
			par.set(y, x);
			sz.set(x, sz[x] + sz[y]);
		}
	}
	void undo() {
		if (hist.empty()) return;
		bool united = hist.back();
		if (united) {
			par.undo();
			sz.undo();
			components++;
		}
		hist.pop_back();
	}
	int cc() { return components; }
};

template<class T> class unionfind_with_potential {
	vector<int> par;
	vector<int> sz;
	vector<vector<int>> group;
public:
	vector<T> pot;
	unionfind_with_potential(int n) {
		par = vector<int>(n);
		for (int i = 0; i < n; i++) par[i] = i;
		sz = vector<int>(n, 1);
		pot = vector<T>(n);
		group.resize(n);
		rep(i, n) group[i].push_back(i);
	}
	int find(int x) {
		if (par[x] == x) return x;
		else return par[x] = find(par[x]);
	}
	int size(int x) { return sz[find(x)]; }
	bool same(int x, int y) { return find(x) == find(y); }
	bool unite(int xv, int yv, T d) {
		int x = find(xv);
		int y = find(yv);
		if (x == y) return pot[xv] == pot[yv] + d;
		T c = pot[xv] - pot[yv] - d;
		if (sz[x] < sz[y]) {
			par[x] = y;
			sz[y] += sz[x];
			for (int e : group[x]) {
				group[y].push_back(e);
				pot[e] -= c;
			}
			group[x].clear();
		}
		else {
			par[y] = x;
			sz[x] += sz[y];
			for (int e : group[y]) {
				group[x].push_back(e);
				pot[e] += c;
			}
			group[y].clear();
		}
		return true;
	}
};

// ABC235Ex, ARC098F, CGR30-E
// (u, v) を結ぶ重み w の辺
// 重み昇順にマージ, 同じ重みは同時に
// 元の頂点 i -> 頂点 i に
// 親は子より index 大
// 子があるなら 2 個以上あります (重み distinct ならちょうど 2 個)
// (森構造, マージされたときの重み, 根の集合) を返す
tuple<vector<vector<int>>, vector<int>, vector<int>> UF_merge(int n, vector<array<int, 3>> E) {
	vector<int> ind(n);
	rep(i, n) ind[i] = i;
	sort(all(E), [&](array<int, 3> s, array<int, 3> t) { return s[2] < t[2]; });
	vector<vector<int>> G(n);
	vector<int> W(n);
	unionfind uf(n);
	rep(i, ssize(E)) {
		vector<P<int>> p;
		{
			auto [u, v, w] = E[i];
			if (!uf.same(u, v)) {
				p.push_back({ uf.find(u), uf.find(v) });
				uf.unite(u, v);
			}
		}
		while (i + 1 < ssize(E) && E[i][2] == E[i + 1][2]) {
			i++;
			auto [u, v, w] = E[i];
			if (!uf.same(u, v)) {
				p.push_back({ uf.find(u), uf.find(v) });
				uf.unite(u, v);
			}
		}
		if (p.empty()) continue;
		vector<int> use;
		for (auto& [u, v] : p) {
			use.push_back(u);
			use.push_back(v);
		}
		sort(all(use));
		use.erase(unique(all(use)), use.end());
		int use_sz = ssize(use);
		unionfind tmp_uf(use_sz);
		for (auto& [u, v] : p) {
			u = lower_bound(all(use), u) - use.begin();
			v = lower_bound(all(use), v) - use.begin();
			tmp_uf.unite(u, v);
		}
		vector tmp_vec(use_sz, vector<int>());
		rep(j, use_sz) tmp_vec[tmp_uf.find(j)].push_back(j);
		int cur = ssize(G);
		rep(j, use_sz) {
			if (tmp_uf.find(j) == j) {
				vector<int> new_G;
				for (int k : tmp_vec[j]) {
					int v = use[k];
					G[ind[v]].push_back(cur);
					new_G.push_back(ind[v]);
				}
				G.push_back(new_G);
				W.push_back(E[i][2]);
				ind[uf.find(use[j])] = cur;
				cur++;
			}
		}
	}
	vector<int> root;
	rep(i, n) if (uf.find(i) == i) root.push_back(ind[i]);
	return { G, W, root };
}