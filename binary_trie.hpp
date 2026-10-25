template<class S, const int D, bool UNIQUE> class binary_trie {
	static_assert(is_integral_v<S>);
	static_assert(0 < D && D < numeric_limits<S>::digits);
	struct Node { int ch[2]{}; int c = 0; };
	vector<Node> data = vector<Node>(2); // 0 : null, 1 : root
	S ofs = 0;
public:
	void reserve(int n) { data.reserve(n); }
	int size() { return data[1].c; }
	bool empty() { return data[1].c == 0; }
	void insert(S x) {
		assert(0 <= x && x < (S(1) << D));
		x ^= ofs;
		array<int, D + 1> path;
		path[D] = 1;
		int cur = 1;
		for (int i = D - 1; i >= 0; i--) {
			int ind = x >> i & 1;
			if (!data[cur].ch[ind]) {
				int nxt = ssize(data);
				data.emplace_back();
				data[cur].ch[ind] = nxt;
			}
			cur = data[cur].ch[ind];
			path[i] = cur;
		}
		if (UNIQUE && data[cur].c) return;
		rep(i, D + 1) data[path[i]].c++;
	}
	void erase(S x) {
		assert(0 <= x && x < (S(1) << D));
		x ^= ofs;
		array<int, D + 1> path;
		path[D] = 1;
		int cur = 1;
		for (int i = D - 1; i >= 0; i--) {
			int ind = x >> i & 1;
			int nxt = data[cur].ch[ind];
			if (data[nxt].c == 0) return;
			cur = nxt;
			path[i] = cur;
		}
		rep(i, D + 1) data[path[i]].c--;
	}
	S min_xor(S x) {
		assert(0 <= x && x < (S(1) << D));
		x ^= ofs;
		int cur = 1;
		assert(data[cur].c != 0);
		S res = 0;
		for (int i = D - 1; i >= 0; i--) {
			int ind = x >> i & 1;
			int nxt = data[cur].ch[ind];
			if (data[nxt].c) cur = nxt;
			else {
				cur = data[cur].ch[ind ^ 1];
				res += S(1) << i;
			}
		}
		return res;
	}
	S max_xor(S x) {
		assert(0 <= x && x < (S(1) << D));
		x ^= ofs;
		int cur = 1;
		assert(data[cur].c != 0);
		S res = 0;
		for (int i = D - 1; i >= 0; i--) {
			int ind = (x >> i & 1) ^ 1;
			int nxt = data[cur].ch[ind];
			if (data[nxt].c) {
				cur = nxt;
				res += S(1) << i;
			}
			else cur = data[cur].ch[ind ^ 1];
		}
		return res;
	}
	void xor_all(S x) { 
		assert(0 <= x && x < (S(1) << D)); 
		ofs ^= x;
	}
	// x の個数
	int count(S x) {
		assert(0 <= x && x < (S(1) << D));
		x ^= ofs;
		int cur = 1;
		for (int i = D - 1; i >= 0; i--) {
			int ind = x >> i & 1;
			cur = data[cur].ch[ind];
			if (cur == 0) return 0;
		}
		return data[cur].c;
	}
	// [0, x) の個数
	int count_range(S x) {
		if (x <= 0) return 0;
		if (x >= (S(1) << D)) return size();
		int cur = 1;
		int res = 0;
		for (int i = D - 1; i >= 0; i--) {
			int u = ofs >> i & 1;
			if (x >> i & 1) {
				res += data[data[cur].ch[u]].c;
				cur = data[cur].ch[u ^ 1];
				if (cur == 0) break;
			}
			else {
				cur = data[cur].ch[u];
				if (cur == 0) break;
			}
		}
		return res;
	}
	int count_range(S l, S r) {
		assert(0 <= l && l <= r && r <= (S(1) << D));
		return count_range(r) - count_range(l);
	}
	// 昇順で k 番目
	S kth_number(int k) {
		assert(0 <= k && k < size());
		int cur = 1;
		S res = 0;
		for (int i = D - 1; i >= 0; i--) {
			int ind = ofs >> i & 1;
			if (data[data[cur].ch[ind]].c <= k) {
				k -= data[data[cur].ch[ind]].c;
				cur = data[cur].ch[ind ^ 1];
				res += S(1) << i;
			}
			else cur = data[cur].ch[ind];
		}
		return res;
	}
	S lower_bound(S x) {
		int k = count_range(x);
		if (k == size()) return S(1) << D;
		return kth_number(k);
	}
	S upper_bound(S x) { return lower_bound(x + 1); }
};