template<class T> int lis(vector<T>& a) {
	int n = ssize(a);
	auto MX = numeric_limits<T>::max();
	for (auto& e : a) assert(e < MX);
	vector<T> tmp(n, MX);
	rep(i, n) {
		int idx = lower_bound(all(tmp), a[i]) - tmp.begin();
		tmp[idx] = a[i];
	}
	int res = lower_bound(all(tmp), MX) - tmp.begin();
	return res;
}

template<class T> int lds(vector<T>& a) {
	int n = ssize(a);
	auto MN = numeric_limits<T>::min();
	for (auto& e : a) assert(e > MN);
	vector<T> tmp(n, MN);
	rep(i, n) {
		int idx = lower_bound(all(tmp), a[i], greater<>()) - tmp.begin();
		tmp[idx] = a[i];
	}
	int res = lower_bound(all(tmp), MN, greater<>()) - tmp.begin();
	return res;
}
