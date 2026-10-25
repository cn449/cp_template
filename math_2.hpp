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

/*
ll extgcd(ll a, ll b, ll& p, ll& q) {
	if (b == 0) { p = 1; q = 0; return a; }
	ll d = extgcd(b, a % b, q, p);
	q -= a / b * p;
	return d;
}

ll my_lcm(ll x, ll y) {
	ll g = gcd(x, y);
	return mul(x, y / g);
}

P<ll> crt(P<ll> f, P<ll> g) {
	ll p, q;
	ll d = extgcd(f.first, g.first, p, q);
	assert((g.second - f.second) % d == 0);
	ll md = my_lcm(f.first, g.first);
	ll tmp = (g.second - f.second) / d * p % (g.first / d);
	ll r = (f.second + f.first * tmp) % md;
	if (r < 0) r += md;
	return { md, r };
}
*/

/*vector<int> primitive_root(A);

struct calc_primitive_root {
	calc_primitive_root() {
		for (int p = 2; p < A; p++) {
			if (pf[p] != p) continue;
			ll tmp = p - 1;
			vector<int> prime_factor;
			while (tmp > 1) {
				ll q = pf[tmp];
				prime_factor.push_back(q);
				while (tmp % q == 0) tmp /= q;
			}
			while (true) {
				ll u = rng(1, p - 1);
				bool ok = true;
				for (ll q : prime_factor) {
					if (modpow(u, (p - 1) / q, p) == 1) {
						ok = false;
						break;
					}
				}
				if (ok) {
					primitive_root[p] = u;
					break;
				}
			}
		}
	}
}calc_primitive_root_;*/