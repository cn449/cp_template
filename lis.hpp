int lis(vector<int> &a) {
	int n = ssize(a);
	vector<int> tmp(n, inf);
	rep(i, n) {
		int idx = lower_bound(all(tmp), a[i]) - tmp.begin();
		tmp[idx] = a[i];
	}
	int res = lower_bound(all(tmp), inf) - tmp.begin();
	return res;
}