using sps = vector<mint>;
// using Poly = vector<mint>;

// 3^N を落とすために N <= 20 などで出題されがちなので TL / ML に注意 
// N = 20, mint で 1700ms 程度
// N = 20 で int / mint だと 300 MB 弱、long long だと 600MB 弱
// a = b のときは片方のみゼータ変換する / a, b は長さ popcount + 1 までしか持たないなどの工夫もできる

void sps_zeta(vector<Poly>& a) {
	int n = countr_zero<uint>(ssize(a));
	rep(i, n) {
		rep(j, 1 << n) {
			if (j >> i & 1) rep(k, n + 1) a[j][k] += a[j - (1 << i)][k];
			else j += (1 << i) - 1;
		}
	}
}

void sps_mobius(vector<Poly>& a) {
	int n = countr_zero<uint>(ssize(a));
	rep(i, n) {
		rep(j, 1 << n) {
			if (j >> i & 1) rep(k, n + 1) a[j][k] -= a[j - (1 << i)][k];
			else j += (1 << i) - 1;
		}
	}
}

sps subset_convolution(sps& a, sps& b) {
	int n = countr_zero<uint>(ssize(a));
	vector<Poly> f(1 << n), g(1 << n);
	rep(i, 1 << n) {
		f[i].resize(n + 1);
		f[i][popcount<uint>(i)] = a[i];
		g[i].resize(n + 1);
		g[i][popcount<uint>(i)] = b[i];
	}
	sps_zeta(f); sps_zeta(g);
	rep(i, 1 << n) {
		for (int j = n; j >= 0; j--) {
			f[i][j] *= g[i][0];
			for (int k = 1; k <= j; k++) f[i][j] += f[i][j - k] * g[i][k];
		}
		g[i].clear(); g[i].shrink_to_fit();
	}
	sps_mobius(f);
	vector<mint> res(1 << n);
	rep(i, 1 << n) res[i] = f[i][popcount<uint>(i)];
	return res;
}

sps sps_exp(sps& a) {
	int n = countr_zero<uint>(ssize(a));
	vector f(1 << n, vector<mint>());
	rep(i, 1 << n) {
		f[i].resize(n + 1);
		f[i][popcount<uint>(i)] = a[i];
	}
	sps_zeta(f);
	rep(i, 1 << n) f[i] = sparse_exp(f[i]);
	sps_mobius(f);
	vector<mint> res(1 << n);
	rep(i, 1 << n) res[i] = f[i][popcount<uint>(i)];
	return res;
}

sps sps_log(sps& a) {
	int n = countr_zero<uint>(ssize(a));
	vector f(1 << n, vector<mint>());
	rep(i, 1 << n) {
		f[i].resize(n + 1);
		f[i][popcount<uint>(i)] = a[i];
	}
	sps_zeta(f);
	rep(i, 1 << n) f[i] = sparse_log(f[i]);
	sps_mobius(f);
	vector<mint> res(1 << n);
	rep(i, 1 << n) res[i] = f[i][popcount<uint>(i)];
	return res;
}

sps composite_poly_sps(Poly p, sps a) {
	int m = ssize(p), n = countr_zero<uint>ssize(a));
	mint const_term = a[0];
	a[0] = 0;
	vector<mint> pw(m + 1);
	pw[0] = 1;
	rep(i, m) pw[i + 1] = pw[i] * const_term;
	vector<mint> g(n + 1);
	rep(i, m) rep(j, n + 1) if (i >= j) g[j] += p[i] * pw[i - j] * comb(i, j);

	vector f(1 << n, vector<mint>());
	rep(i, 1 << n) {
		f[i].resize(n + 1);
		f[i][popcount<uint>(i)] = a[i];
	}
	sps_zeta(f);
	rep(i, 1 << n) {
		f[i] = composition_slow(g, f[i]);
		f[i].resize(n + 1);
	}
	sps_mobius(f);
	sps res(1 << n);
	rep(i, 1 << n) res[i] = f[i][popcount<uint>(i)];
	return res;
}