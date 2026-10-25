ll binary_gcd(ll a, ll b) {
	a = abs(a), b = abs(b);
	if (a == 0) return b;
	if (b == 0) return a;
	int c = countr_zero<ull>(a | b);
	a >>= countr_zero<ull>(a);
	do {
		b >>= countr_zero<ull>(b);
		if (a > b) swap(a, b);
		b -= a;
	} while (b);
	return a << c;
}