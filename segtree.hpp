template<class T, auto op, auto e> class segtree {
	int n;
	vector<T> data;
public:
	segtree(int _n) {
		n = 1;
		while (n < _n + 2) n <<= 1;
		data = vector<T>(2 * n, e());
	}
	segtree(vector<T> vec) {
		int _n = vec.size();
		n = 1;
		while (n < _n + 2) n <<= 1;
		data = vector<T>(2 * n, e());
		for (int i = 0; i < _n; i++) data[i + n] = vec[i];
		for (int i = n - 1; i >= 1; i--) data[i] = op(data[i << 1], data[i << 1 | 1]);
	}
	void change(int i, T x) {
		i += n;
		data[i] = x;
		while (i > 1) {
			i >>= 1;
			data[i] = op(data[i << 1], data[i << 1 | 1]);
		}
	}
	void add(int i, T x) { change(i, data[i + n] + x); }
	T get(int a, int b) {
		T left = e(); T right = e();
		a += n; b += n;
		while (a < b) {
			if (a & 1) left = op(left, data[a++]);
			if (b & 1) right = op(data[--b], right);
			a >>= 1; b >>= 1;
		}
		return op(left, right);
	}
	T get_all() { return data[1]; }
	T operator[](int i) { return data[i + n]; }
};