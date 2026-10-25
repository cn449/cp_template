template<class S, const int D> class wavelet_matrix {
	int n;
	vector<S> a;
	array<vector<int>, D> data;
	int count_0(int d, int l, int r) { return r - l - data[d][r] + data[d][l]; }
	int count_1(int d, int l, int r) { return data[d][r] - data[d][l]; }
	void succ_0(int d, int& l, int& r) { l = count_0(d, 0, l), r = count_0(d, 0, r); }
	void succ_1(int d, int& l, int& r) { l = count_0(d, 0, n) + count_1(d, 0, l), r = count_0(d, 0, n) + count_1(d, 0, r); }
public:
	wavelet_matrix(const vector<S>& _a) : n(ssize(_a)), a(_a) {
		vector<S> na(n);
		for (int d = D - 1; d >= 0; d--) {
			data[d].resize(n + 1);
			rep(i, n) data[d][i + 1] = data[d][i] + (a[i] >> d & 1);
			int ind = 0;
			rep(rp, 2) rep(i, n) if ((a[i] >> d & 1) == rp) na[ind++] = a[i];
			swap(a, na);
		}
	}
	S kth_number(int l, int r, int k) {
		assert(0 <= l && l < r && r <= n);
		assert(0 <= k && k < r - l);
		S res = 0;
		for (int d = D - 1; d >= 0; d--) {
			int c0 = count_0(d, l, r);
			if (k < c0) succ_0(d, l, r);
			else {
				k -= c0;
				res += S(1) << d;
				succ_1(d, l, r);
			}
		}
		return res;
	}
	int range_freq(int l, int r, S upper) {
		assert(0 <= l && l <= r && r <= n);
		if (upper >= S(1) << D) return r - l;
		int res = 0;
		for (int d = D - 1; d >= 0; d--) {
			if (upper >> d & 1) {
				res += count_0(d, l, r);
				succ_1(d, l, r);
			}
			else succ_0(d, l, r);
		}
		return res;
	}
	int range_freq(int l, int r, S lower, S upper) {
		assert(0 <= l && l <= r && r <= n);
		if (lower > upper) return 0;
		return range_freq(l, r, upper) - range_freq(l, r, lower);
	}
};

template<class S, class T, const int D> class wavelet_matrix_sum {
	int n;
	vector<S> a;
	vector<T> b;
	array<vector<int>, D> data;
	array<vector<T>, D> sb;
	int count_0(int d, int l, int r) { return r - l - data[d][r] + data[d][l]; }
	int count_1(int d, int l, int r) { return data[d][r] - data[d][l]; }
	void succ_0(int d, int& l, int& r) { l = count_0(d, 0, l), r = count_0(d, 0, r); }
	void succ_1(int d, int& l, int& r) { l = count_0(d, 0, n) + count_1(d, 0, l), r = count_0(d, 0, n) + count_1(d, 0, r); }
public:
	wavelet_matrix_sum(const vector<S>& _a, const vector<T> &_b) : n(ssize(_a)), a(_a), b(_b) {
		vector<S> na(n);
		vector<T> nb(n);
		for (int d = D - 1; d >= 0; d--) {
			data[d].resize(n + 1);
			sb[d].resize(n + 1);
			rep(i, n) data[d][i + 1] = data[d][i] + (a[i] >> d & 1);
			int ind = 0;
			rep(rp, 2) rep(i, n) if ((a[i] >> d & 1) == rp) {
				na[ind] = a[i];
				nb[ind] = b[i];
				ind++;
			}
			swap(a, na);
			swap(b, nb);
			rep(i, n) sb[d][i + 1] = sb[d][i] + b[i];
		}
	}
	T sum(int l, int r, S upper) {
		T res = 0;
		for (int d = D - 1; d >= 0; d--) {
			if (upper >> d & 1) {
				int l0 = count_0(d, 0, l), r0 = count_0(d, 0, r);
				res += sb[d][r0] - sb[d][l0];
				succ_1(d, l, r);
			}
			else succ_0(d, l, r);
		}
		return res;
	}
	T sum(int l, int r, S lower, S upper) {
		return sum(l, r, upper) - sum(l, r, lower);
	}
};

template<class S, class T, const int D> class wavelet_bit {
	int n;
	vector<S> a;
	array<vector<int>, D> data;
	array<fenwick_tree<T>, D> ds;
	int count_0(int d, int l, int r) { return r - l - data[d][r] + data[d][l]; }
	int count_1(int d, int l, int r) { return data[d][r] - data[d][l]; }
	void succ_0(int d, int& l, int& r) { l = count_0(d, 0, l), r = count_0(d, 0, r); }
	void succ_1(int d, int& l, int& r) { l = n - count_1(d, l, n), r = n - count_1(d, r, n); }
public:
	wavelet_bit(const vector<S>& _a) : n(ssize(_a)), a(_a) {
		vector<S> na(n);
		for (int d = D - 1; d >= 0; d--) {
			data[d].resize(n + 1);
			ds[d] = fenwick_tree<T>(n);
			rep(i, n) data[d][i + 1] = data[d][i] + (a[i] >> d & 1);
			int ind = 0;
			rep(rp, 2) rep(i, n) if ((a[i] >> d & 1) == rp) na[ind++] = a[i];
			swap(a, na);
		}
	}
	void add(int pos, T val) {
		for (int d = D - 1; d >= 0; d--) {
			if (data[d][pos] == data[d][pos + 1]) pos = count_0(d, 0, pos);
			else pos = n - count_1(d, pos, n);
			ds[d].add(pos, val);
		}
	}
	T sum(int l, int r, S upper) {
		T res = 0;
		for (int d = D - 1; d >= 0; d--) {
			if (upper >> d & 1) {
				int l0 = count_0(d, 0, l), r0 = count_0(d, 0, r);
				res += ds[d].sum(l0, r0);
				succ_1(d, l, r);
			}
			else succ_0(d, l, r);
		}
		return res;
	}
};