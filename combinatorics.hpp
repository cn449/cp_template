// CF1105-D
// [L, R] から C 個選んだ積の和 (同じ数 ok, 順序なし)
// L = 1 のときは第 2 種スターリング数 S(R + C, R)
mint sum_of_prod(int L, int R, int C) {
	mint res = 0;
	int D = R - L;
	rep(i, D + 1) res += mint(D % 2 == i % 2 ? 1 : -1) * mint(L + i).pow(C + D) * facinv[i] * facinv[D - i];
	return res;
}

// ARC202-D
// + に n 回, - に m 回
// prefix sum がすべて [L, R] の範囲
// N = n + m, D = R - L として O(N/D)
// 上下限制約の一方しかないとき O(1)
mint generalized_dyck_path(int n, int m, int L, int R) {
	if (!(L <= 0 && 0 <= R)) return 0;
	if (!(L <= n - m && n - m <= R)) return 0;
	mint res = 0;
	int t = R - L + 2, w = n - R - 1;
	for (int i = 0; n + i * t <= n + m; i++) res += comb(n + m, n + i * t);
	for (int i = -1; n + i * t >= 0; i--) res += comb(n + m, n + i * t);
	for (int i = 0; w + i * t <= n + m; i++) res -= comb(n + m, w + i * t);
	for (int i = -1; w + i * t >= 0; i--) res -= comb(n + m, w + i * t);
	return res;
}

vector<int> pf; // 最小素因数, pf[0] = pf[1] = 0
vector<int> prime; // 素数列

// 複数回やると壊れる（それはそう）
void sieve(int A) {
	pf.resize(A);
	for (int i = 2; i < A; i++) {
		if (pf[i] != 0) continue;
		// prime.push_back(i);
		for (int j = i; j < A; j += i) if (pf[j] == 0) pf[j] = i;
	}
}

// 0^d, 1^d,..., {n-1}^d 列挙
vector<mint> enum_power(int n, int d) {
	assert(n <= ssize(pf));
	vector<mint> res(n);
	if (n >= 1) res[0] = d == 0 ? 1 : 0;
	if (n >= 2) res[1] = 1;
	for (int i = 2; i < n; i++) {
		if (pf[i] == i) res[i] = mint(i).pow(d);
		else {
			int p = pf[i], q = i / p;
			res[i] = res[p] * res[q];
		}
	}
	return res;
}

// O(n)
// 嘘で、線形篩にしてないから O(n loglog n) だ
mint Bernoulli(int n) {
	assert(n >= 0);
	if (n == 0) return 1;
	if (n == 1) {
		assert(mod % 2 == 1);
		return (mod - 1) / 2;
	}
	if (n % 2 == 1) return 0;
	mint f0 = 0;
	rep(i, n + 1) f0 += inv[i + 1];
	vector<mint> A(n);
	rep(i, n) A[i] = comb(n, i + 1) * (i % 2 == 0 ? 1 : -1);
	A = naive_conv(A, fps{ 1, -1 });
	A[0] += 1;
	auto pw = enum_power(n + 1, n);
	mint p = f0;
	mint res = p * pw[0];
	rep(k, n) {
		mint t = inv[k + 1] * A[k];
		p = p - t;
		res += p * pw[k + 1];
	}
	return res;
}