ll safe_sqrt(ll n) {
	assert(n >= 0);
	if (n == 0) return 0;
	ll w = sqrtl(n);
	while (w * w > n) w--;
	while ((w + 1) * (w + 1) <= n) w++;
	return w;
}

// floor(N / x) として現れる値の管理
// N = 10
// quo : [10,5,3,2,1]
// itv : [[1,2),[2,3),[3,4),[4,6),[6,11)] -> [1,2,3,4,6,11]
// access_0 : [1,2,3,4,5,6,7,8,9,10] -> [4,3,2,1,1,0,0,0,0,0]
// access_1 : [1,2,3,4,5,6,7,8,9,10] -> [0,1,2,3,3,4,4,4,4,4]

class quotients {
	ll N;
	vector<ll> quo;
	vector<ll> itv;
	int large, small;
public:
	quotients() {}
	quotients(ll N_) : N(N_) {
		assert(N <= ten(16));
		ll L = 1;
		while (L <= N) {
			ll R = N / (N / L) + 1;
			quo.push_back(N / L);
			itv.push_back(L);
			L = R;
		}
		itv.push_back(L);
		small = safe_sqrt(N);
		large = ssize(quo) - small;
	}
	ll get(int i) { return quo[i]; }
	P<ll> get_itv(int i) { return { itv[i], itv[i + 1] }; }
	int size() { return ssize(quo); }
	int access_0(ll x) {
		if (x <= small) return size() - x;
		else return N / x - 1;
	}
	int access_1(ll x) { return access_0(N / x); }
};