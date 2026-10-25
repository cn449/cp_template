template<class T, auto op, auto e> class Sparsetable {
	int n;
	vector<vector<T>> data;
public:
	Sparsetable() {}
	Sparsetable(vector<T> v) {
		n = ssize(v);
		data.push_back(v);
		for (int i = 2; i <= n; i <<= 1) {
			vector<T> w(n - i + 1);
			for (int j = 0; j <= n - i; j++) {
				w[j] = op(data.back()[j], data.back()[j + i / 2]);
			}
			data.push_back(w);
		}
	}
	T prod(int l, int r) {
		if (l >= r) return e();
		int u = bit_floor<uint>(r - l);
		int d = countr_zero<uint>(u);
		return op(data[d][l], data[d][r - u]);
	}
};