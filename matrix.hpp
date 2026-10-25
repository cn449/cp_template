template<class T> vector<vector<T>> mat_mul(vector<vector<T>>& a, vector<vector<T>>& b) {
	int n = ssize(a), m = n == 0 ? 0 : ssize(a[0]);
	int m2 = ssize(b), l = m2 == 0 ? 0 : ssize(b[0]);
	assert(m == m2);
	rep(i, n) assert(ssize(a[i]) == m);
	rep(i, m) assert(ssize(b[i]) == l);
	vector res(n, vector<T>(l));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			for (int k = 0; k < l; k++) {
				res[i][k] += a[i][j] * b[j][k];
			}
		}
	}
	return res;
}

template<class T> vector<vector<T>> mat_pow(vector<vector<T>> a, ll b) {
	int n = ssize(a);
	rep(i, n) assert(ssize(a[i]) == n);
	vector res(n, vector<T>(n));
	for (int i = 0; i < n; i++) res[i][i] = 1;
	while (b) {
		if (b & 1) res = mat_mul(res, a);
		a = mat_mul(a, a);
		b >>= 1;
	}
	return res;
}

template<class T> int mat_rank(vector<vector<T>> a) {
	int n = ssize(a), m = n == 0 ? 0 : ssize(a[0]);
	rep(i, n) assert(ssize(a[i]) == m);
	if (m > n) {
		vector b(m, vector<T>(n));
		rep(i, n) rep(j, m) b[j][i] = a[i][j];
		a = move(b); swap(n, m);
	}
	int res = 0;
	rep(i, n) {
		int f = -1;
		for (int j = res; j < m; j++) if (a[i][j] != 0) { f = j; break; }
		if (f == -1) continue;
		if (res != f) rep(j, n) swap(a[j][res], a[j][f]);
		T co = a[i][res]; T coinv = (T)1 / co;
		rep(j, m) a[i][j] *= coinv;
		for (int j = i + 1; j < n; j++) {
			T r = a[j][res];
			rep(k, m) a[j][k] -= a[i][k] * r;
		}
		res++;
	}
	return res;
}

template<class T> T det(vector<vector<T>> a) {
	int n = ssize(a);
	rep(i, n) assert(ssize(a) == n);
	T res = 1;
	rep(i, n) {
		int f = -1;
		for (int j = i; j < n; j++) if (a[j][i] != 0) { f = j; break; }
		if (f == -1) return 0;
		if (i != f) { swap(a[i], a[f]); res *= -1; }
		T co = a[i][i]; T coinv = (T)1 / co; res *= co;
		rep(j, n) a[i][j] *= coinv;
		for (int j = i + 1; j < n; j++) {
			T r = a[j][i];
			rep(k, n) a[j][k] -= a[i][k] * r;
		}
	}
	return res;
}

template<class T> pair<bool, vector<vector<T>>> mat_inv(vector<vector<T>> a) {
	int n = ssize(a);
	rep(i, n) assert(ssize(a[i]) == n);
	vector res(n, vector<T>(n));
	rep(i, n) res[i][i] = 1;
	rep(i, n) {
		int f = -1;
		for (int j = i; j < n; j++) if (a[j][i] != 0) { f = j; break; }
		if (f == -1) return { false, {} };
		if (i != f) { swap(a[i], a[f]); swap(res[i], res[f]); }
		T co = a[i][i]; T coinv = (T)1 / co;
		rep(j, n) { a[i][j] *= coinv; res[i][j] *= coinv; }
		rep(j, n) {
			if (j == i) continue;
			T r = a[j][i];
			rep(k, n) { a[j][k] -= a[i][k] * r; res[j][k] -= res[i][k] * r; }
		}
	}
	return { true, res };
}

// O(nm^2) 時間
// ax = v の解 x を返す
// first : 解空間の次元（存在しないとき -1）
// second : 0 番目の要素は解の 1 つ、1 番目以降は解空間の基底
// 行とか列が 0 になるときに気を付ける（別途処理が必要がち？）
template<class T> pair<int, vector<vector<T>>> linear_equation(vector<vector<T>> a, vector<T> v, bool unique = false) {
	int n = ssize(a), m = n == 0 ? 0 : ssize(a[0]);
	assert(ssize(v) == n);
	int cur = 0;
	rep(j, m) {
		int f = -1;
		for (int i = cur; i < n; i++) if (a[i][j] != 0) { f = i; break; }
		if (f == -1) continue;
		if (cur != f) { swap(a[cur], a[f]); swap(v[cur], v[f]); }
		T co = a[cur][j]; T coinv = (T)1 / co;
		rep(k, m) a[cur][k] *= coinv;
		v[cur] *= coinv;
		rep(i, n) {
			if (i == cur) continue;
			T r = a[i][j];
			rep(k, m) a[i][k] -= a[cur][k] * r;
			v[i] -= v[cur] * r;
		}
		cur++;
		if (cur == n) break;
	}
	vector<T> res(m);
	vector<int> idx(m, -1);
	rep(i, n) {
		int f = -1;
		rep(j, m) if (a[i][j] != 0) { f = j; break; }
		if (f == -1 && v[i] != 0) { return { -1,{} }; }
		if (f != -1) { res[f] = v[i]; idx[f] = i; }
	}
	if (unique) return { 0, {res} };
	vector<vector<T>> vec = { res };
	int dim = 0;
	rep(j, m) {
		if (idx[j] != -1) continue;
		vector<T> w(m); w[j] = -1;
		rep(i, m) if (idx[i] != -1) w[i] = a[idx[i]][j];
		dim++; vec.push_back(w);
	}
	return { dim,vec };
}

template<class T> vector<T> characteristic_poly(vector<vector<T>> a) {
	int n = ssize(a);
	rep(i, n - 1) {
		int f = -1;
		for (int j = i + 1; j < n; j++) {
			if (a[j][i] != 0) {
				f = j;
				break;
			}
		}
		if (f != -1) {
			rep(j, n) swap(a[i + 1][j], a[f][j]);
			rep(j, n) swap(a[j][f], a[j][i + 1]);
			for (int j = i + 2; j < n; j++) {
				T coeff = a[j][i] / a[i + 1][i];
				rep(k, n) a[j][k] -= a[i + 1][k] * coeff;
				rep(k, n) a[k][i + 1] += a[k][j] * coeff;
			}
		}
	}
	vector dp(n + 1, vector<T>(n + 1));
	dp[0][0] = 1;
	for (int i = 1; i <= n; i++) {
		For(j, 1, i + 1) dp[i][j] = dp[i - 1][j - 1];
		rep(j, i) dp[i][j] -= a[i - 1][i - 1] * dp[i - 1][j];
		T coeff = 1;
		for (int j = i - 2; j >= 0; j--) {
			coeff *= a[j + 1][j];
			rep(k, j + 1) dp[i][k] -= a[j][i - 1] * coeff * dp[j][k];
		}
	}
	return dp[n];
}

mint count_undirected_spanning_trees(int n, vector<vector<int>> E) {
	assert(ssize(E) == n);
	rep(i, n) assert(ssize(E[i]) == n);
	vector M(n - 1, vector<mint>(n - 1));
	rep(i, n - 1) rep(j, n) M[i][i] += E[i][j];
	rep(i, n - 1) rep(j, n - 1) M[i][j] -= E[i][j];
	mint res = det(M);
	return res;
}

// 根に向かう方向
mint count_directed_spanning_trees(int n, vector<vector<int>> E, int root) {
	assert(ssize(E) == n);
	rep(i, n) assert(ssize(E[i]) == n);
	vector M(n, vector<mint>(n));
	rep(i, n) rep(j, n) M[i][i] += E[i][j];
	rep(i, n) rep(j, n) M[i][j] -= E[i][j];
	vector A(n - 1, vector<mint>(n - 1));
	rep(i, n - 1) rep(j, n - 1) A[i][j] = M[i < root ? i : i + 1][j < root ? j : j + 1];
	mint res = det(A);
	return res;
}

mint count_eulerian_circuits(int n, vector<vector<int>> E) {
	assert(ssize(E) == n);
	rep(i, n) assert(ssize(E[i]) == n);
	vector<int> in(n), out(n);
	rep(i, n) rep(j, n) {
		in[j] += E[i][j];
		out[i] += E[i][j];
	}
	rep(i, n) rep(j, n) if (in[i] != out[i]) return 0;
	vector<int> vec;
	rep(i, n) {
		bool use = false;
		rep(j, n) if (E[i][j] > 0) use = true;
		if (use) vec.push_back(i);
	}
	int sz = ssize(vec);
	if (sz == 0) return 1; // 0 辺のときは答え 1 とする
	vector M(sz - 1, vector<mint>(sz - 1));
	rep(i, sz - 1) rep(j, sz - 1) {
		if (i == j) M[i][j] = out[vec[j]];
		M[i][j] -= E[vec[i]][vec[j]];
	}
	mint ans = det(M);
	for (int i : vec) ans *= fac[out[i] - 1];
	return ans;
}

// distinguish : 多重辺を区別しますか？
// 区別しない出題が多い気がするなぁ
mint count_eulerian_path(int n, vector<vector<int>> E, bool distinguish = false) {
	mint sum = 0;
	vector<int> in(n), out(n);
	rep(i, n) rep(j, n) {
		in[j] += E[i][j];
		out[i] += E[i][j];
		sum += E[i][j];
	}
	int s = -1, t = -1;
	rep(i, n) {
		int d = in[i] - out[i];
		if (abs(d) > 1) return 0;
		if (d == 1) {
			if (t == -1) t = i;
			else return 0;
		}
		if (d == -1) {
			if (s == -1) s = i;
			else return 0;
		}
	}
	mint coeff = 1;
	if (!distinguish) rep(i, n) rep(j, n) coeff *= facinv[E[i][j]];
	if (s == -1 && t == -1) coeff *= sum;
	else E[t][s]++;
	return count_eulerian_circuits(n, E) * coeff;
}

vector<P<int>> general_matching(int n, vector<P<int>> E) {
	vector a(n, vector<mint>(n));
	for (auto [u, v] : E) {
		assert(u != v);
		mint x = rng(1, mod - 1);
		a[u][v] = x;
		a[v][u] = -x;
	}
	auto memo_a = a;
	vector<int> vertex;
	int rk = 0;
	rep(i, n) {
		int f = -1;
		rep(j, n) if (a[i][j] != 0) { f = j; break; }
		if (f == -1) continue;
		rep(j, n) swap(a[j][rk], a[j][f]);
		assert(a[i][rk] != 0);
		mint coeff = mint(1) / a[i][rk];
		rep(j, n) {
			if (i == j) continue;
			mint c = coeff * a[j][rk];
			rep(k, n) a[j][k] -= c * a[i][k];
		}
		vertex.push_back(i);
		rk++;
	}
	vector A(rk, vector<mint>(rk));
	rep(i, rk) rep(j, rk) A[i][j] = memo_a[vertex[i]][vertex[j]];
	auto [exist, B] = mat_inv(A);
	assert(exist);
	vector<P<int>> res;
	rep(_, rk / 2) {
		int sz = ssize(A);
		int u = -1, v = -1;
		rep(i, sz) {
			rep(j, sz) if (A[i][j] != 0 && B[i][j] != 0) { u = i, v = j; break; }
			if (u != -1) break;
		}
		assert(u != -1 && v != -1);
		res.push_back({ vertex[u], vertex[v] });
		assert(u < v);
		swap(vertex[sz - 1], vertex[v]); swap(vertex[sz - 2], vertex[u]);
		swap(A[sz - 1], A[v]); swap(A[sz - 2], A[u]);
		rep(i, sz) {
			swap(A[i][sz - 1], A[i][v]);
			swap(A[i][sz - 2], A[i][u]);
		}
		swap(B[sz - 1], B[v]); swap(B[sz - 2], B[u]);
		rep(i, sz) {
			swap(B[i][sz - 1], B[i][v]);
			swap(B[i][sz - 2], B[i][u]);
		}
		swap(B[sz - 1][sz - 2], B[sz - 2][sz - 1]);
		B[sz - 1][sz - 2] = mint(1) / B[sz - 1][sz - 2];
		B[sz - 2][sz - 1] = mint(1) / B[sz - 2][sz - 1];
		rep(i, sz - 2) rep(j, sz - 2) {
			B[i][j] -= B[i][sz - 2] * B[sz - 2][sz - 1] * B[sz - 1][j];
			B[i][j] -= B[i][sz - 1] * B[sz - 1][sz - 2] * B[sz - 2][j];
		}
		vertex.resize(sz - 2);
		A.resize(sz - 2); B.resize(sz - 2);
		rep(i, sz - 2) {
			A[i].resize(sz - 2);
			B[i].resize(sz - 2);
		}
	}
	return res;
}