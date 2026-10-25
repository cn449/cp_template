// https://judge.yosupo.jp/submission/269450
// CGR28-F, ABC435-F
// 各頂点の親の index の列を返す
template<class T> vector<int> CartesianTree(vector<T> &a) {
	int n = ssize(a);
	vector<int> par(n, -1);
	vector<pair<int, T>> vec;
	rep(i, n) {
		int pos = -1;
		while (!vec.empty() && vec.back().second > a[i]) {
			pos = vec.back().first;
			vec.pop_back();
		}
		if (pos != -1) par[pos] = i;
		if (!vec.empty()) par[i] = vec.back().first;
		vec.push_back({ i, a[i] });
	}
	return par;
}