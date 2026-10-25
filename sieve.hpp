vector<int> pf; // 最小素因数, pf[0] = pf[1] = 0
vector<int> prime; // 素数列

// 複数回やると壊れる（それはそう）
void sieve(int A) {
	pf.resize(A);
	for (int i = 2; i < A; i++) {
		if (pf[i] != 0) continue;
		prime.push_back(i);
		for (int j = i; j < A; j += i) if (pf[j] == 0) pf[j] = i;
	}
}