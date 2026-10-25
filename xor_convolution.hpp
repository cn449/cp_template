mint inv2 = (mod + 1) / 2;

void fwt(vector<mint>& f) {
	int n = ssize(f);
	for (int i = 1; i < n; i <<= 1) {
		for (int j = 0; j < n; j++) {
			if ((j & i) == 0) {
				mint x = f[j], y = f[j | i];
				f[j] = x + y, f[j | i] = x - y;
			}
		}
	}
}

void ifwt(vector<mint>& f) {
	int n = ssize(f);
	for (int i = 1; i < n; i <<= 1) {
		for (int j = 0; j < n; j++) {
			if ((j & i) == 0) {
				mint x = f[j], y = f[j | i];
				f[j] = (x + y) * inv2, f[j | i] = (x - y) * inv2;
			}
		}
	}
}