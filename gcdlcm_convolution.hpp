vector<int> pf; // 最小素因数, pf[0] = pf[1] = 0
vector<int> prime; // 素数列

void sieve(int A) {
	pf.resize(A);
	for (int i = 2; i < A; i++) {
		if (pf[i] != 0) continue;
		prime.push_back(i);
		for (int j = i; j < A; j += i) if (pf[j] == 0) pf[j] = i;
	}
}

template<class T> void gcd_zeta(vector<T>& a) {
	int n = ssize(a);
	for (int p : prime) {
		for (int i = (n - 1) / p; i >= 1; i--) {
			a[i] += a[i * p];
		}
	}
}

template<class T> void gcd_mobius(vector<T>& a) {
	int n = ssize(a);
	for (int p : prime) {
		for (int i = 1; i <= (n - 1) / p; i++) {
			a[i] -= a[i * p];
		}
	}
}

template<class T> vector<T> gcd_convolution(vector<T> a, vector<T> b) {
	int n = max(ssize(a), ssize(b));
	a.resize(n); b.resize(n);
	gcd_zeta(a); gcd_zeta(b);
	rep(i, n) a[i] *= b[i];
	gcd_mobius(a);
	return a;
}

template<class T> void lcm_zeta(vector<T>& a) {
	int n = ssize(a);
	for (int p : prime) {
		for (int i = 1; i <= (n - 1) / p; i++) {
			a[i * p] += a[i];
		}
	}
}

template<class T> void lcm_mobius(vector<T>& a) {
	int n = ssize(a);
	for (int p : prime) {
		for (int i = (n - 1) / p; i >= 1; i--) {
			a[i * p] -= a[i];
		}
	}
}

// lcm を取ると大きくなることがある → sz まで拡張
template<class T> vector<T> lcm_convolution(vector<T> a, vector<T> b, int sz = -1) {
	int n = max({ int(ssize(a)), int(ssize(b)), sz });
	a.resize(n); b.resize(n);
	lcm_zeta(a); lcm_zeta(b);
	rep(i, n) a[i] *= b[i];
	lcm_mobius(a);
	return a;
}