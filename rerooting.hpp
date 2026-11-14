template<class S, class T = S>
struct Rerooting {
    int n, m = 0;
    vector<vector<pair<int, int>>> G;
    Rerooting(int n) : n(n), G(n) {}
    Rerooting(const vector<vector<int>>& G_) : Rerooting(ssize(G_)) {
        rep(i, n) for (int j : G_[i]) if (i < j) add_edge(i, j);
    }
    void add_edge(int u, int v) {
        assert(0 <= u && u < n && 0 <= v && v < n);
        G[u].push_back({v, m});
        G[v].push_back({u, m});
        m++;
    }
    template<class Op, class E, class Vertex, class Edge>
    vector<T> solve(Op op, E e, Vertex vertex, Edge edge, int root = 0) const {
        assert(n > 0);
        assert(0 <= root && root < n);
        assert(m == n - 1);
        S id = e();
        vector<int> par(n, -1), pe(n, -1), order;
        order.reserve(n);
        par[root] = root;
        order.push_back(root);
        rep(i, ssize(order)) {
            int v = order[i];
            for (auto [to, eid] : G[v]) if (par[to] == -1) {
                par[to] = v;
                pe[to] = eid;
                order.push_back(to);
            }
        }
        assert(ssize(order) == n);
        vector<S> down(n, id), up(n, id);
        for (int i = n - 1; i >= 0; i--) {
            int v = order[i];
            S x = id;
            for (auto [to, eid] : G[v]) if (to != par[v]) x = op(x, down[to]);
            if (v != root) down[v] = edge(vertex(x, v), v, par[v], pe[v]);
        }
        vector<T> ans(n);
        vector<S> pref;
        for (int v : order) {
            pref.clear();
            pref.push_back(id);
            for (auto [to, eid] : G[v]) pref.push_back(op(pref.back(), to == par[v] ? up[v] : down[to]));
            ans[v] = vertex(pref.back(), v);
            S suf = id;
            for (int i = ssize(G[v]) - 1; i >= 0; i--) {
                auto [to, eid] = G[v][i];
                if (to != par[v]) up[to] = edge(vertex(op(pref[i], suf), v), v, to, eid);
                suf = op(to == par[v] ? up[v] : down[to], suf);
            }
        }
        return ans;
    }
};
