namespace rdv {
	random_device seed_gen;
	mt19937_64 engine(seed_gen());
	ll rnum(ll r) { return engine() % r; } // [0, r)
	ll rnum(ll l, ll r) { return rnum(r - l) + l; } // [l, r)
	ll rng(ll l, ll r) { return rnum(l, r + 1); } // [l, r]
	double rng01() { return engine() * pow(2, -64); }
	template<class T> void shuf(vector<T>& v) { shuffle(all(v), engine); }
	void shuf(string& s) { shuffle(all(s), engine); }
	vector<vector<int>> make_vec(int n, int l, int r) {
		if (n == 0) return { {} };
		vector<vector<int>> pre = make_vec(n - 1, l, r);
		vector<vector<int>> res;
		for (vector<int> v : pre) {
			for (int i = l; i < r; i++) {
				vector<int> w = v;
				w.push_back(i);
				res.push_back(w);
			}
		}
		return res;
	}
	vector<vector<int>> make_vec(int n, int r) { return make_vec(n, 0, r); }
	vector<vector<int>> make_vec_inc(int n, int l, int r) {
		if (n == 0) return { {} };
		vector<vector<int>> pre = make_vec_inc(n - 1, l, r);
		vector<vector<int>> res;
		for (vector<int> v : pre) {
			int b = l - 1;
			if (!v.empty()) b = v.back();
			for (int i = b + 1; i < r; i++) {
				vector<int> w = v;
				w.push_back(i);
				res.push_back(w);
			}
		}
		return res;
	}
	vector<vector<int>> make_vec_inc(int n, int r) { return make_vec_inc(n, 0, r); }
	vector<vector<int>> make_vec_sum(int n, int sum) {
		if (n == 0) return { {} };
		vector<vector<int>> pre = make_vec_sum(n - 1, sum);
		vector<vector<int>> res;
		for (vector<int> v : pre) {
			int lim = sum;
			for (int w : v) lim -= w;
			for (int i = 0; i <= lim; i++) {
				vector<int> w = v;
				w.push_back(i);
				res.push_back(w);
			}
		}
		return res;
	}
	vector<vector<int>> make_rand_tree(int n) {
		vector<vector<int>> g(n, vector<int>());
		for (int i = 1; i < n; i++) {
			int v = rng(0, i - 1);
			g[i].push_back(v);
			g[v].push_back(i);
		}
		return g;
	}
}