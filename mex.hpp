template<class T> int mex(vector<T> &a) {
	int n = ssize(a);
	vector<bool> used(n + 1);
	rep(i, n) if (0 <= a[i] && a[i] <= n) used[a[i]] = true;
	rep(i, n + 1) if (!used[i]) return i;
}