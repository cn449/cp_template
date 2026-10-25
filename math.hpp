ll modinv(ll a, ll m) {
	if (m == 1) return a;
	a %= m;
	if (a == 0) abort();
	ll b = m, u = 1, v = 0;
	while (b) {
		ll t = a / b;
		a -= t * b; swap(a, b);
		u -= t * v; swap(u, v);
	}
	u %= m;
	if (u < 0) u += m;
	return u;
}

ll modpow(ll a, ll b, ll m) {
	ll ans = 1;
	a %= m;
	while (b) {
		if (b & 1) ans = ans * a % m;
		a = a * a % m;
		b >>= 1;
	}
	return ans;
}

// mod m1 で r1, mod m2 で r2
P<ll> crt(ll m1, ll r1, ll m2, ll r2) {
	r1 %= m1; r2 %= m2;
	if (m1 == -1 or m2 == -1) return { -1, -1 };
	if (m1 == m2) {
		if (r1 == r2) return { m1,r1 };
		else return { -1,-1 };
	}
	ll g = gcd(m1, m2);
	if (r1 % g != r2 % g) return { -1, -1 };
	ll res = r1 % g;
	r1 /= g; m1 /= g;
	r2 /= g; m2 /= g;
	if (m2 == 1) {
		swap(r1, r2);
		swap(m1, m2);
	}
	ll v = r2 - r1;
	if (v < 0) {
		ll d = -v;
		v += (d + m2 - 1) / m2 * m2;
	}
	v *= modinv(m1, m2);
	v %= m2;
	ll m = m1 * m2 * g;
	return { m, ((r1 + v * m1) * g + res) % m };
}

// mod first で second
P<ll> crt_many(vector<P<ll>> f) {
	ll m = 1, r = 0;
	for (auto p : f) {
		P<ll> nxt = crt(m, r, p.first, p.second);
		m = nxt.first, r = nxt.second;
	}
	return { m,r };
}

// a * b mod m
// a, b, m は 64 bit に収まる
ll mul_long(ll a, ll b, ll m) {
	return lint(a) * b % m;
}

// a ^ b mod m
// a, b, m は 64 bit に収まる
ll pow_long(ll a, ll b, ll m) {
	ll res = 1;
	while (b) {
		if (b & 1) res = mul_long(res, a, m);
		a = mul_long(a, a, m);
		b >>= 1;
	}
	return res;
}

// https://judge.yosupo.jp/submission/271045
// 定数倍微妙？ 340ms で 10^18 オーダーを 10^5 個
bool is_prime(ll n) {
	if (n <= 1) return false;
	if (n == 2 or n == 3 or n == 5 or n == 7) return true;
	if (n % 2 == 0 or n % 3 == 0 or n % 5 == 0 or n % 7 == 0) return false;
	vector<ll> a;
	if (n < (1ull << 32)) a = { 2,7,61 };
	else a = { 2,325,9375,28178,450775,9780504,1795265022 };
	ll d = n - 1, s = 0;
	while (d % 2 == 0) { d /= 2; s++; }
	for (ll b : a) {
		if (n <= b) return true;
		bool flag = false;
		ll x = pow_long(b, d, n);
		if (x == 1 or x == n - 1) flag = true;
		else {
			rep(_, s) {
				x = mul_long(x, x, n);
				if (x == n - 1) flag = true;
			}
		}
		if (!flag) return false;
	}
	return true;
}

// 合成数に対し非自明な因数を返す
ll find_factor(ll n) {
	if (n % 2 == 0) return 2;
	if (n % 3 == 0) return 3;
	if (n % 5 == 0) return 5;
	if (n % 7 == 0) return 7;
	for (ll c = 1;; c++) {
		auto f = [&](ll x) -> ll { return ((lint)x * x + c) % n; };
		ll x = 0, y = 0;
		while (true) {
			x = f(x);
			y = f(f(y));
			ll g = gcd(x - y, n);
			if (g == n) break;
			else if (g != 1) return g;
		}
	}
}

// 素因数分解
vector<ll> factorize(ll n) {
	if (n == 1) return {};
	if (is_prime(n)) return { n };
	vector<ll> res;
	vector<ll> vec = { n };
	while (!vec.empty()) {
		ll v = vec.back(); vec.pop_back();
		ll d0 = find_factor(v);
		ll d1 = v / d0;
		if (is_prime(d0)) res.push_back(d0);
		else if (d0 != 1) vec.push_back(d0);
		if (is_prime(d1)) res.push_back(d1);
		else if (d1 != 1) vec.push_back(d1);
	}
	sort(all(res));
	return res;
}

// 原始根を 1 つ求める
ll primitive_root(ll p) {
	if (p == 2) return 1;
	if (p == 3) return 2;
	if (p == 998244353) return 3;
	vector<ll> f = factorize(p - 1);
	f.erase(unique(all(f)), f.end());
	while (true) {
		ll v = rng(2, p - 1);
		bool valid = true;
		for (auto q : f) {
			if (pow_long(v, (p - 1) / q, p) == 1) valid = false;
		}
		if (valid) return v;
	}
}

// ax = b (mod m) を x について解く、解なしなら -1
ll lin_eq(ll a, ll b, ll m) {
	a %= m; b %= m;
	if (a == 0 && b == 0) return 0;
	ll g = gcd(a, m);
	if (b % g != 0) return -1;
	a /= g; b /= g; m /= g;
	ll x = b * modinv(a, m) % m;
	return x;
}

// n : mod p で平方元ですか？
bool is_square_modp(ll n, ll p) {
	n %= p;
	if (n < 0) n += p;
	if (n == 0 or n == 1) return true;
	return modpow(n, (p - 1) / 2, p) == 1;
}

// u : 非平方元として適当に選んでくる
ll mod_sqrt_(ll n, ll p, ll u) {
	n %= p;
	if (n < 0) n += p;
	if (n == 0 or n == 1) return n;
	int f = 0;
	ll q = p - 1;
	while (q % 2 == 0) {
		q /= 2;
		if (modpow(n, q, p) != 1) break;
		f++;
	}
	if (((p - 1) / (1 << f)) % 2 == 1) { return modpow(n, (q + 1) / 2, p); }
	ll v = modpow(u, 1 << f, p);
	ll res = modpow(u, 1 << (f - 1), p) * mod_sqrt_(n * modinv(v, p) % p, p, u) % p;
	return res;
}

// x^2 = n mod p なる x を返す, なければ -1
ll mod_sqrt(ll n, ll p) {
	n %= p;
	if (n < 0) n += p;
	if (n == 0 or n == 1) return n;
	if (!is_square_modp(n, p)) return -1;
	ll u = 1;
	while (is_square_modp(u, p)) u = rng(1, p - 1);
	return mod_sqrt_(n, p, u);
}