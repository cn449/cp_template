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

using fps = vector<mint>;
using Poly = vector<mint>;

void dft(Poly& f) { internal::butterfly(f); }
void idft(Poly& f) {
	internal::butterfly_inv(f);
	mint coeff = mint(1) / ssize(f);
	for (auto& e : f) e *= coeff;
}

// 愚直掛け算
fps naive_conv(fps f, fps g) {
	int n = ssize(f), m = ssize(g);
	fps res(n + m - 1);
	rep(i, n) rep(j, m) res[i + j] += f[i] * g[j];
	return res;
}

fps operator +(fps a, fps b) {
	int n = ssize(a), m = ssize(b);
	fps res(max(n, m));
	rep(i, n) res[i] += a[i];
	rep(i, m) res[i] += b[i];
	return res;
}

fps operator -(fps a, fps b) {
	int n = ssize(a), m = ssize(b);
	fps res(max(n, m));
	rep(i, n) res[i] += a[i];
	rep(i, m) res[i] -= b[i];
	return res;
}

fps operator *(fps a, fps b) {
	if (mod == mod998) return convolution(a, b);
	return naive_conv(a, b);
}

fps convolution_large(const fps& a, const fps& b) {
	int n = ssize(a), m = ssize(b);
	const int B = 1 << 22;
	int sa = (n + B - 1) / B, sb = (m + B - 1) / B;
	vector<fps> f(sa, fps(2 * B)), g(sb, fps(2 * B));
	rep(i, n) f[i >> 22][i & (B - 1)] = a[i];
	rep(i, m) g[i >> 22][i & (B - 1)] = b[i];
	rep(i, sa) dft(f[i]);
	rep(i, sb) dft(g[i]);
	vector<fps> h(sa + sb - 1, fps(2 * B));
	rep(i, sa) rep(j, sb) rep(k, 2 * B) h[i + j][k] += f[i][k] * g[j][k];
	rep(i, sa + sb - 1) idft(h[i]);
	fps res(n + m - 1);
	rep(i, sa + sb - 1) rep(j, 2 * B) if (h[i][j] != 0) res[(i << 22) + j] += h[i][j];
	return res;
}

// f の逆元を mod x^deg で求める
// 何も指定しないと f の与えた次数まで
fps inverse(fps f, int deg = -1) {
	assert(f[0] != 0);
	if (deg != -1) f.resize(deg);
	int n = ssize(f);
	fps res(n);
	res[0] = mint(1) / mint(f[0]);
	fps g, h;
	for (int k = 1; k < n; k <<= 1) {
		g.resize(4 * k); h.resize(4 * k);
		rep(i, 4 * k) g[i] = 0, h[i] = 0;
		rep(i, min(n, 2 * k)) { g[i] = res[i]; h[i] = f[i]; }
		dft(g); dft(h);
		rep(i, 4 * k) { g[i] *= g[i]; g[i] *= h[i]; }
		idft(g);
		for (int i = k; i < n && i < 2 * k; i++) res[i] -= g[i];
	}
	res.resize(n);
	return res;
}

fps operator /(fps f, fps g) {
	int d = max(ssize(f), ssize(g));
	auto res = f * inverse(g, d);
	res.resize(d);
	return res;
}

void operator +=(fps& a, fps b) { a = a + b; }
void operator -=(fps& a, fps b) { a = a - b; }
void operator *=(fps& a, fps b) { a = a * b; }
void operator /=(fps& a, fps b) { a = a / b; }

// FPS としての割り算（愚直）
// g が sparse のときにも使える
fps naive_div(fps f, fps g) {
	assert(!g.empty()); assert(g[0] != 0);
	int n = ssize(f), m = ssize(g);
	mint c = mint(1) / mint(g[0]);
	for (auto& e : f) e *= c;
	for (auto& e : g) e *= c;
	vector<int> idx;
	for (int i = 1; i < m; i++) if (g[i] != 0) idx.push_back(i);
	vector<mint> res(n);
	rep(i, n) {
		res[i] = f[i];
		for (int j : idx) {
			if (i + j >= ssize(f)) break;
			f[i + j] -= res[i] * g[j];
		}
	}
	return res;
}

// 微分
fps diff(fps f) {
	int n = ssize(f);
	fps res(n - 1);
	rep(i, n - 1) res[i] = f[i + 1] * (i + 1);
	return res;
}

// 積分
fps inte(fps f) {
	int n = ssize(f);
	fps res(n + 1);
	rep(i, n) res[i + 1] = f[i] * inv[i + 1];
	return res;
}

fps log(fps f) {
	assert(f[0] == 1);
	if (ssize(f) == 1) return { 0 };
	fps res = diff(f) / f;
	res.insert(res.begin(), 0);
	res.resize(ssize(f));
	For(i, 1, ssize(f)) res[i] *= inv[i];
	return res;
}

fps exp(fps f) {
	assert(f[0] == 0);
	int n = ssize(f);
	fps res = { 1, 0 };
	fps g = { 1 }, h = { 1 };
	for (int k = 1; k < n; k <<= 1) {
		h.resize(2 * k);
		swap(g, res);
		For(i, k, min(n, 2 * k)) h[i] = f[i];
		res = g * (h - log(g));
		res.resize(4 * k);
	}
	res.resize(n);
	return res;
}

fps pow(fps f, ll a) {
	int n = ssize(f);
	if (a == 0) { fps O(n); O[0] = 1; return O; }
	int s = -1;
	rep(i, n) { if (f[i] != 0) { s = i; break; } }
	if (s == -1) return { fps(n) };
	mint fs = f[s];
	mint fsinv = mint(1) / mint(fs);
	for (mint& v : f) v *= fsinv;
	fps fnew(ssize(f) - s);
	for (int i = s; i < ssize(f); i++) fnew[i - s] = f[i];
	fps g = log(fnew);
	for (mint& v : g) v *= a;
	fps r = exp(g);
	mint w = fs.pow(a);
	for (mint& v : r) v *= w;
	r.resize(n);
	fps res(n);
	ll L = safe_mul(s, a);
	for (ll i = L; i < n; i++) res[i] = r[i - s * a];
	return res;
}

pair<bool, fps> sqrt(fps f) {
	int n = ssize(f);
	int s = -1;
	rep(i, n) { if (f[i] != 0) { s = i; break; } }
	if (s == -1) return { true, fps(n) };
	if (s % 2 == 1) return { false, {} };
	mint fs = f[s];
	for (mint& v : f) v /= fs;
	fps fnew;
	for (int i = s; i < ssize(f); i++) fnew.push_back(f[i]);
	fps g = log(fnew);
	for (mint& v : g) v *= inv[2];
	fps r = exp(g);
	ll w = mod_sqrt(fs.val(), mod);
	if (w == -1) return { false, {} };
	for (mint& v : r) v *= w;
	r.resize(n);
	fps res(n);
	rep(i, n) if (i >= s / 2) res[i] = r[i - s / 2];
	return { true, res };
}

fps sparse_log(fps f) {
	auto g = diff(f);
	auto h = naive_div(g, f);
	return inte(h);
}

fps sparse_exp(fps f) {
	assert(f[0] == 0);
	int n = ssize(f);
	fps g(n); // df/dx
	rep(i, n - 1) g[i] = f[i + 1] * (i + 1);
	vector<int> idx;
	rep(i, n) if (g[i] != 0) idx.push_back(i);
	fps res(n);
	res[0] = 1;
	rep(i, n - 1) {
		mint sum = 0;
		for (int j : idx) {
			if (i < j) break;
			sum += g[j] * res[i - j];
		}
		res[i + 1] = sum * inv[i + 1];
	}
	return res;
}

fps sparse_pow(fps f, ll a) {
	int n = ssize(f);
	if (a == 0) { fps O(n); O[0] = 1; return O; }
	ll s = 0;
	rep(i, n) { if (f[i] == 0) s++; else break; }
	ll ofs = safe_mul(s, a);
	if (ofs >= n) return fps(n);
	ll d = n - ofs;
	fps fnew(n - s);
	rep(i, n - s) fnew[i] = f[i + s];
	n -= s;
	f = move(fnew);
	vector<int> idx_f;
	rep(i, n) if (f[i] != 0) idx_f.push_back(i);
	fps g(n - 1);
	rep(i, n - 1) g[i] = f[i + 1] * (i + 1);
	vector<int> idx_g;
	rep(i, n - 1) if (g[i] != 0) idx_g.push_back(i);
	fps r(d + ofs);
	r[ofs] = f[0].pow(a);
	fps h(d), t(d);
	mint f0inv = mint(1) / mint(f[0]);
	rep(i, d - 1) {
		for (int j : idx_g) {
			if (i < j) break;
			h[i] += g[j] * r[ofs + i - j] * a;
		}
		for (int j : idx_f) {
			if (j == 0) continue;
			if (i < j) break;
			h[i] -= f[j] * t[i - j];
		}
		t[i] = h[i] * f0inv;
		r[ofs + i + 1] = t[i] * inv[i + 1];
	}
	return r;
}

pair<bool, fps> sparse_sqrt(fps f) {
	int n = ssize(f);
	ll s = 0;
	rep(i, n) { if (f[i] == 0) s++; else break; }
	if (s == n) return { true, fps(n) };
	if (s % 2 == 1) return { false, {} };
	ll ofs = s / 2;
	ll d = n - ofs;
	fps fnew(n - s);
	rep(i, n - s) fnew[i] = f[i + s];
	n -= s;
	f = move(fnew);
	ll sqf0 = mod_sqrt(f[0].val(), mod);
	if (sqf0 == -1) return { false, {} };
	vector<int> idx_f;
	rep(i, n) if (f[i] != 0) idx_f.push_back(i);
	fps g(n - 1);
	rep(i, n - 1) g[i] = f[i + 1] * (i + 1);
	vector<int> idx_g;
	rep(i, n - 1) if (g[i] != 0) idx_g.push_back(i);
	fps r(d + ofs);
	r[ofs] = sqf0;
	fps h(d), t(d);
	mint f0inv = mint(1) / mint(f[0]);
	rep(i, d - 1) {
		for (int j : idx_g) {
			if (i < j) break;
			h[i] += g[j] * r[ofs + i - j] * inv[2];
		}
		for (int j : idx_f) {
			if (j == 0) continue;
			if (i < j) break;
			h[i] -= f[j] * t[i - j];
		}
		t[i] = h[i] * f0inv;
		r[ofs + i + 1] = t[i] * inv[i + 1];
	}
	return { true, r };
}

// 多項式としての割り算
// (商, 余り) の組を返す
P<Poly> div_poly(Poly f, Poly g) {
	while (!f.empty() && f.back() == 0) f.pop_back();
	while (!g.empty() && g.back() == 0) g.pop_back();
	int n = ssize(f), m = ssize(g);
	if (n < m) { return { {0}, f }; }
	Poly revf = f, revg = g;
	reverse(all(revf)); reverse(all(revg));
	Poly q = revf / revg;
	q.resize(n - m + 1);
	reverse(all(q));
	Poly r = f - g * q;
	r.resize(m - 1);
	return { q, r };
}

// ABC245D
// 多項式としての割り算(愚直)
// (商, 余り) の組を返す
P<Poly> naive_div_poly(Poly f, Poly g) {
	while (!g.empty() && g.back() == 0) g.pop_back();
	assert(!g.empty());
	int n = ssize(f), m = ssize(g);
	n--; m--; // deg f = n, deg g = m
	mint inv = mint(1) / mint(g.back());
	int sz = n - m;
	vector<mint> q(sz + 1);
	for (int i = n; i >= m; i--) {
		q[i - m] = f[i] * inv;
		rep(j, m + 1) f[i - m + j] -= q[i - m] * g[j];
	}
	vector<mint> r(m);
	rep(i, m) r[i] = f[i];
	return { q, r };
}

// ABC315Ex
struct RelaxedConvolution {
	int M, q;
	vector<mint> a, b, c;
	vector<vector<vector<mint>>> f, g;
	RelaxedConvolution(int _n) {
		M = 0, q = 0;
		int n = 1;
		while (n < _n + 5) { n <<= 1; M++; }
		a.resize(1 << M); b.resize(1 << M); c.resize(1 << M); f.resize(M); g.resize(M);
	}
	mint get(mint x, mint y) {
		assert(q <= (1 << M) - 2);
		a[q] = x, b[q] = y;
		for (int i = 0; i < M; i++) {
			int p = q + 2, d = 1 << i;
			if (p == 2 * d) {
				vector<mint> ax(d), bx(d);
				for (int j = 0; j < d; j++) ax[j] = a[j + d - 1], bx[j] = b[j + d - 1];
				vector<mint> r = ax * bx;
				f[i].push_back(r); g[i].push_back(r);
			}
			else if (p > 2 * d && p % d == 0) {
				vector<mint> ax(d), bx(d), ay(d), by(d);
				for (int j = 0; j < d; j++) ax[j] = a[j + p - d - 1], bx[j] = b[j + d - 1], ay[j] = a[j + d - 1], by[j] = b[j + p - d - 1];
				f[i].push_back(ax * bx); g[i].push_back(ay * by);
			}
		}
		mint res = 0;
		for (int i = 0; i < M; i++) {
			int d = 1 << i, s = 2 * d - 2, t = 2 * s;
			if (q < s) continue;
			int idx = (q - s) / d;
			res += f[i][idx][q - s - idx * d];
			if (idx != 0) res += g[i][idx][q - s - idx * d];
			idx--;
			if (idx >= 0 && q <= t + idx * d) {
				res += f[i][idx][q - s - idx * d];
				if (idx != 0) res += g[i][idx][q - s - idx * d];
			}
		}
		q++;
		return res;
	}
};

// CF783E
// c[i - j] += a[i] * b[j]
// a と同じ長さが返ってくる
fps rev_conv(fps a, fps b) {
	int m = ssize(b);
	reverse(all(b));
	fps c = a * b;
	int sz = ssize(c);
	fps res(sz - m + 1);
	rep(i, sz - m + 1) res[i] = c[i + m - 1];
	return res;
}

// f(x + c) を多項式として計算
Poly taylor_shift(Poly f, mint c) {
	int n = ssize(f);
	rep(i, n) f[i] *= fac[i];
	Poly g(n);
	mint coeff = 1;
	rep(i, n) {
		g[i] = facinv[i] * coeff;
		coeff *= c;
	}
	auto h = rev_conv(f, g);
	rep(i, n) h[i] *= facinv[i];
	return h;
}

// 多項式 f の 1 点 v での値を求める
mint Eval(Poly& f, mint v) {
	mint res = 0, d = 1;
	for (mint e : f) {
		res += d * e;
		d *= v;
	}
	return res;
}

// 多項式 f の a[i] たちでの値を求める
vector<mint> MultipointEvaluation(Poly f, vector<mint> a) {
	int m = ssize(a);
	int s = 1;
	while (s < m) s <<= 1;
	vector<fps> v(2 * s);
	rep(i, m) v[s + i] = { -a[i],1 };
	for (int i = m; i < s; i++) v[s + i] = { 1 };
	for (int i = s - 1; i >= 1; i--) v[i] = v[i << 1] * v[i << 1 | 1];
	vector<fps> w(2 * s);
	w[1] = div_poly(f, v[1]).second;
	for (int i = 2; i < s + m; i++) {
		w[i] = div_poly(w[i >> 1], v[i]).second;
		w[i].resize(ssize(v[i]) - 1);
	}
	vector<mint> res(m);
	rep(i, m) res[i] = w[i + s][0];
	return res;
}

// 多項式 f の a,ar,ar^2,...,ar^(m-1) での値を求める
vector<mint> MultipointEvaluation_geometric(Poly f, mint a, mint r, int m) {
	int n = ssize(f);
	if (r == 0) {
		vector<mint> res(m);
		res[0] = Eval(f, a);
		For(i, 1, m) res[i] = f[0];
		return res;
	}
	vector<mint> F(n), G(n + m);
	mint ir = mint(1) / mint(r);
	vector<mint> pwr(n + m + 1), pwir(n + m + 1);
	pwr[0] = 1, pwir[0] = 1;
	rep(i, n + m) {
		pwr[i + 1] = pwr[i] * r;
		pwir[i + 1] = pwir[i] * ir;
	}
	mint coeff = 1;
	rep(i, n) {
		F[i] = f[i] * coeff;
		coeff *= a;
		coeff *= pwir[i];
	}
	coeff = 1;
	rep(i, n + m) {
		G[i] = coeff;
		coeff *= pwr[i];
	}
	auto H = rev_conv(G, F);
	H.resize(m);
	coeff = 1;
	rep(i, m) {
		H[i] *= coeff;
		coeff *= pwir[i];
	}
	return H;
}

fps all_prod(vector<fps>& a) {
	int n = ssize(a);
	if (n == 0) return { 1 };
	auto rec = [&](auto rec, int l, int r) -> fps {
		if (l + 1 == r) return a[l];
		int m = midpoint(l, r);
		return rec(rec, l, m) * rec(rec, m, r);
	};
	return rec(rec, 0, n);
}

P<fps> all_sum_rational(vector<P<fps>>& a) {
	int n = ssize(a);
	if (n == 0) return { {0}, {1} };
	auto rec = [&](auto rec, int l, int r) -> P<fps> {
		if (l + 1 == r) return a[l];
		int m = midpoint(l, r);
		auto [f0, g0] = rec(rec, l, m);
		auto [f1, g1] = rec(rec, m, r);
		return { f0 * g1 + g0 * f1, g0 * g1 };
	};
	return rec(rec, 0, n);
}

fps all_sum_rational(vector<P<fps>>& a, int sz) {
	auto [f, g] = all_sum_rational(a);
	f.resize(sz); g.resize(sz);
	return f / g;
}

// n 次多項式, 0,1,...,n での値が分かっている
// t での値を返す
// ARC033-D
mint interpolate(int n, vector<mint> a, ll t) {
	assert(ssize(a) >= n + 1);
	vector<mint> p(n + 1);
	rep(i, n + 1) {
		p[i] = a[i] * facinv[i] * facinv[n - i];
		if ((n - i) % 2 == 1) p[i] *= -1;
	}
	mint res = 0;
	vector<mint> l(n + 2), r(n + 2);
	l[0] = r[0] = 1;
	rep(i, n + 1) l[i + 1] = l[i] * (t - i);
	rep(i, n + 1) r[i + 1] = r[i] * (t - n + i);
	reverse(all(r));
	rep(i, n + 1) res += l[i] * r[i + 1] * p[i];
	return res;
}

// n 次多項式, 0,1,...,n での値が分かっている
// 係数列を返す, NTT を使わない
// CF901-E
vector<mint> interpolate_coeff_naive(int n, vector<mint> a) {
	assert(ssize(a) >= n + 1);
	// Σ(p[i] * (x - 0)(x - 1)...(x - n) / (x - i))
	vector<mint> p(n + 1);
	rep(i, n + 1) {
		p[i] = a[i] * facinv[i] * facinv[n - i];
		if ((n - i) % 2 == 1) p[i] *= -1;
	}
	vector<mint> res(n + 1);
	vector<mint> prod = { 1 };
	rep(i, n + 1) {
		vector<mint> f = { -i,1 };
		prod = naive_conv(prod, f);
	}
	rep(i, n + 1) {
		vector<mint> f = { -i,1 };
		auto [g, _] = naive_div_poly(prod, f);
		rep(j, n + 1) res[j] += g[j] * p[i];
	}
	return res;
}

// n 次以下で n+1 点与えられる
// f(p[i]) = a[i]
Poly interpolate_coeff(int n, vector<mint> p, vector<mint> a) {
	assert(ssize(a) == n + 1);
	assert(ssize(p) == n + 1);
	vector<Poly> f(n + 1);
	rep(i, n + 1) f[i] = { -p[i], 1 };
	Poly F = all_prod(f);
	F = diff(F);
	vector<mint> denom = MultipointEvaluation(F, p);
	auto rec = [&](auto rec, int L, int R) -> P<Poly> {
		if (L + 1 == R) {
			mint t = a[L] / denom[L];
			return { {t}, {-p[L], 1} };
		}
		int M = midpoint(L, R);
		auto [f0, g0] = rec(rec, L, M);
		auto [f1, g1] = rec(rec, M, R);
		return { f0 * g1 + f1 * g0, g0 * g1 };
	};
	return rec(rec, 0, n + 1).first;
}

// n 次以下で n+1 点与えられる
// f(sr^i) = a[i]
// 点は distinct
// サボりで log 2 個付けている, 定数倍は軽い
Poly interpolate_coeff_geometric(int n, mint s, mint r, vector<mint> a) {
	if (n < 0) return { 0 };
	assert(ssize(a) == n + 1);
	vector<mint> p(n + 1);
	p[0] = s;
	rep(i, n) p[i + 1] = p[i] * r;
	vector<Poly> f(n + 1);
	rep(i, n + 1) f[i] = { -p[i], 1 };
	Poly F = all_prod(f);
	F = diff(F);
	vector<mint> denom = MultipointEvaluation_geometric(F, s, r, n + 1);
	auto rec = [&](auto rec, int L, int R) -> P<Poly> {
		if (L + 1 == R) {
			mint t = a[L] / denom[L];
			return { {t}, {-p[L], 1} };
		}
		int M = midpoint(L, R);
		auto [f0, g0] = rec(rec, L, M);
		auto [f1, g1] = rec(rec, M, R);
		return { f0 * g1 + f1 * g0, g0 * g1 };
	};
	return rec(rec, 0, n + 1).first;
}

vector<mint> berlekamp_massey(const vector<mint>& s) {
	int n = ssize(s);
	vector<mint> C = { 1 }, B = { 1 };
	int L = 0, m = -1;
	mint b = 1;
	rep(i, n) {
		mint d = 0;
		rep(j, L + 1) d += C[j] * s[i - j];
		if (d == 0) continue;
		auto T = C;
		mint coeff = d / b;
		int shift = i - m;
		if (ssize(C) < ssize(B) + shift) C.resize(ssize(B) + shift);
		rep(j, ssize(B)) C[j + shift] -= coeff * B[j];
		if (2 * L <= i) L = i + 1 - L, B = T, b = d, m = i;
	}
	C.erase(C.begin());
	for (auto& e : C) e = -e;
	if (ssize(C) * 2 >= ssize(s)) { cerr << "WARNING : The recurrence is may be incorrect.\n"; }
	return C;
}

// [x^n](f/g)
mint bostan_mori_998(fps f, fps g, ll n) {
	assert(mod == mod998);
	assert(g[0] != 0);
	mint inv2 = (mod + 1) / 2;
	while (n > 0) {
		if (n == 1) {
			f.resize(2);
			g.resize(2);
			return (f[1] * g[0] - f[0] * g[1]) / (g[0] * g[0]);
		}
		if (n & 1) {
			n++;
			f.insert(f.begin(), 0);
		}
		if (ssize(f) > n + 1) f.resize(n + 1);
		if (ssize(g) > n + 1) g.resize(n + 1);
		fps nf, ng;
		uint s0 = ssize(f) + ssize(g) - 1, s1 = ssize(g) + ssize(g) - 1;
		int sz = bit_ceil(max(s0, s1));
		chmax(sz, 2);
		f.resize(sz); g.resize(sz); nf.resize(sz / 2); ng.resize(sz / 2);
		dft(f); dft(g);
		rep(i, sz) f[i] *= g[i ^ 1];
		rep(i, sz / 2) nf[i] = (f[i << 1] + f[i << 1 | 1]) * inv2;
		rep(i, sz / 2) ng[i] = g[i << 1] * g[i << 1 | 1];
		idft(nf); idft(ng);
		n >>= 1;
		f = move(nf), g = move(ng);
		f.resize((s0 + 1) / 2); g.resize((s1 + 1) / 2);
	}
	return f[0] / g[0];
}

mint bostan_mori_slow(fps f, fps g, ll n) {
	while (n > 0) {
		fps h = g;
		for (int i = 1; i < ssize(h); i += 2) h[i] *= -1;
		f *= h; g *= h;
		fps nf, ng;
		nf.reserve((ssize(f) + 1) / 2);
		ng.reserve((ssize(g) + 1) / 2);
		for (int i = 0; (i << 1 | (n & 1)) < ssize(f); i++) nf.push_back(f[i << 1 | (n & 1)]);
		for (int i = 0; i < ssize(g); i += 2) ng.push_back(g[i]);
		n >>= 1;
		f = move(nf), g = move(ng);
	}
	return f[0] / g[0];
}

mint bostan_mori(fps f, fps g, ll n) { return mod == mod998 ? bostan_mori_998(f, g, n) : bostan_mori_slow(f, g, n); }

mint find_nth_term(vector<mint>& a, vector<mint>& c, ll n) {
	if (n < ssize(a)) return a[n];
	int d = ssize(c) - 1;
	while (d >= 0 && c[d] == 0) d--;
	if (d == -1) return 0;
	d++;
	assert(ssize(a) >= d); // 初期条件が足りない！
	fps h(d), g(d + 1);
	g[0] = 1;
	rep(i, d) h[i] = a[i];
	rep(i, d) g[i + 1] = -c[i];
	fps f = h * g;
	f.resize(d);
	return bostan_mori(f, g, n);
}

mint BMBM(vector<mint>& a, ll n) {
	if (n < ssize(a)) return a[n];
	auto rec = berlekamp_massey(a);
	return find_nth_term(a, rec, n);
}

// [x^n](f^i*g) を i = 0,1,...,k-1 で列挙
vector<mint> power_projection(int n, int k, fps f, fps g) {
	f.resize(n + 1); g.resize(n + 1);
	mint const_term = f[0];
	f[0] = 0;
	vector<mint> A(k), B(k);
	int N = bit_ceil<uint>(2 * n + 1), M = 4;
	int S = N * M;
	fps C(S), D(S), F, G;
	D[0] = 1;
	rep(i, n + 1) {
		C[i] = g[i];
		D[N + i] = -f[i];
	}
	mint inv2 = (mod + 1) / 2;
	mint w = mint(3).pow(mod - 1 - (mod - 1) / S);
	while (n > 0) {
		dft(C); dft(D);
		F.resize(S / 2); G.resize(S / 2);
		rep(i, S) C[i] *= D[i ^ 1];
		rep(i, S / 2) {
			if (n & 1) F[i] = C[i << 1] - C[i << 1 | 1];
			else F[i] = C[i << 1] + C[i << 1 | 1];
			F[i] *= inv2;
		}
		rep(i, S / 2) G[i] = D[i << 1] * D[i << 1 | 1];
		if (n & 1) {
			mint coeff = 1;
			int j = 0, b = -1;
			for (int i = 1; i < S / 2; i++) {
				coeff *= w;
				b = (S / 2) >> 1;
				while (j & b) { j ^= b; b >>= 1; }
				j |= b; F[j] *= coeff;
			}
		}
		idft(F); idft(G);
		F.resize(S); G.resize(S);
		rep(j, M) For(i, N / 4, N / 2) F[i + j * (N / 2)] = 0, G[i + j * (N / 2)] = 0;
		N >>= 1; M <<= 1; n >>= 1;
		C = move(F), D = move(G);
		F.clear(); G.clear();
	}
	vector<mint> J(M), K(M);
	rep(i, M) J[i] = C[i * N], K[i] = D[i * N];
	auto H = J / K;
	rep(i, min(int(ssize(H)), k)) A[i] = H[i] * facinv[i];
	B[0] = 1;
	rep(i, k - 1) B[i + 1] = B[i] * const_term;
	rep(i, k) B[i] *= facinv[i];
	auto res = A * B;
	res.resize(k);
	rep(i, k) res[i] *= fac[i];
	return res;
}

// [x^i](f^i*g) を i = 0,1,...,k-1 で列挙
vector<mint> power_projection_2(int k, fps f, fps g) {
	assert(f[0] != 0);
	f.resize(k); g.resize(k);
	fps F = inverse(f);
	F.insert(F.begin(), 0);
	fps G = pow(f, k - 1) * g;
	auto res = power_projection(k - 1, k, F, G);
	reverse(all(res));
	return res;
}

// [x^(i+d)](f^i*g) を i = 0,1,...,k-1 で列挙
// 面倒なので f[0] != 0 とか仮定しちゃおう
vector<mint> power_projection_3(int k, fps f, fps g, int d) {
	assert(f[0] != 0);
	f.resize(k + d); g.resize(k + d);
	assert(d > 0); assert(f[0] != 0);
	fps F = inverse(f);
	F.insert(F.begin(), 0);
	fps G = pow(f, k - 1) * g;
	auto res = power_projection(k - 1 + d, k, F, G);
	reverse(all(res));
	return res;
}

fps compositional_inverse(fps f) {
	int n = ssize(f);
	assert(n >= 2); assert(f[0] == 0); assert(f[1] != 0);
	f.erase(f.begin());
	auto res = power_projection_2(n, inverse(f), { 0,1 });
	For(i, 1, n) res[i] *= inv[i];
	return res;
}

// g(f(x))
fps composition(fps g, fps f) {
	int L = min(ssize(g), ssize(f));
	if (L == 0) return {};
	int n = L - 1, n0 = n;

	auto tdft = [&](fps& a) {
		int s = ssize(a);
		idft(a); reverse(a.begin() + 1, a.end());
		for (auto& x : a) x *= s;
	};
	auto tidft = [&](fps& a) {
		int s = ssize(a);
		reverse(a.begin() + 1, a.end()); dft(a);
		mint is = mint(s).inv();
		for (auto& x : a) x *= is;
	};

	f.resize(L);
	g = taylor_shift(g, f[0]);
	f[0] = 0;
	g.resize(L);
	int N = bit_ceil<uint>(2 * n + 1), N0 = N, M = 4;
	int S = N * M;
	fps C, D(S), F, G;
	D[0] = 1;
	rep(i, L) D[N + i] = -f[i];
	mint inv2 = (mod + 1) / 2;
	mint w = mint(3).pow(mod - 1 - (mod - 1) / S);
	vector<fps> hist;
	while (n > 0) {
		dft(D);
		hist.push_back(move(D));
		auto& E = hist.back();
		G.resize(S / 2);
		rep(i, S / 2) G[i] = E[i << 1] * E[i << 1 | 1];
		idft(G);
		G.resize(S);
		rep(j, M) For(i, N / 4, N / 2) G[i + j * (N / 2)] = 0;
		N >>= 1; M <<= 1; n >>= 1;
		D = move(G); G.clear();
	}
	fps K(M), I(M);
	rep(i, M) K[i] = D[i * N];
	I[0] = 1;
	I = I / K;
	I.resize(M);
	F.clear(); F.resize(M);
	rep(i, min(L, M)) F[i] = g[i];
	reverse(all(F));
	G = F * I;
	C.clear(); C.resize(S);
	rep(i, M) C[i * N] = G[M - 1 - i];
	for (int t = ssize(hist) - 1; t >= 0; t--) {
		n = n0 >> t; N = N0 >> t; M = 4 << t;
		auto& E = hist[t];
		F = move(C);
		rep(j, M) For(i, N / 4, N / 2) F[i + j * (N / 2)] = 0;
		F.resize(S / 2);
		tidft(F);
		if (n & 1) {
			mint coeff = 1;
			int j = 0, b = -1;
			for (int i = 1; i < S / 2; i++) {
				coeff *= w;
				b = (S / 2) >> 1;
				while (j & b) { j ^= b; b >>= 1; }
				j |= b; F[j] *= coeff;
			}
		}
		C.resize(S);
		rep(i, S / 2) {
			C[i << 1] = F[i] * inv2;
			C[i << 1 | 1] = F[i] * inv2 * ((n & 1) ? -1 : 1);
		}
		rep(i, S) C[i] *= E[i ^ 1];
		tdft(C);
	}
	fps res(L);
	rep(i, L) res[i] = C[L - 1 - i];
	return res;
}