// CF736-D2
ll binary_gcd(ll a, ll b) {
	a = abs(a), b = abs(b);
	int c = 0;
	while (a > 0 && b > 0) {
		int x = countr_zero<ull>(a), y = countr_zero<ull>(b);
		if (x == 0 && y == 0) {
			if (a < b) swap(a, b);
			a -= b;
		}
		else {
			c += min(x, y);
			a >>= x; b >>= y;
		}
	}
	if (a == 0) return b << c;
	if (b == 0) return a << c;
}
