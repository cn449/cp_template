// CF740-D
// CGR30-F2

#include <optional>

template<class Key, class S, auto op, auto e, class F, auto mapping, auto composition, auto id, auto shift_mapping> struct lazy_shift_index_treap {
	struct Node {
		Key key;
		S val, prod;
		F lz;
		bool has_lz;
		Key key_lz;
		bool has_key_lz;
		uint64_t pri;
		int sz;
		Node* l, * r;

		Node(const Key& k, const S& v, uint64_t p)
			: key(k), val(v), prod(v), lz(id()), has_lz(false),
			key_lz(Key()), has_key_lz(false),
			pri(p), sz(1), l(nullptr), r(nullptr) {}
	};

	Node* root = nullptr;

	static uint64_t rng() {
		static mt19937_64 mt((uint64_t)chrono::steady_clock::now().time_since_epoch().count());
		return mt();
	}

	static int node_size(Node* t) {
		return t ? t->sz : 0;
	}

	static S prod_all(Node* t) {
		return t ? t->prod : e();
	}

	static void pull(Node* t) {
		if (!t) return;
		t->sz = 1 + node_size(t->l) + node_size(t->r);
		t->prod = op(prod_all(t->l), op(t->val, prod_all(t->r)));
	}

	static void apply_val(Node* t, const F& f) {
		if (!t) return;
		t->val = mapping(f, t->val);
		t->prod = mapping(f, t->prod);
		if (t->has_lz) t->lz = composition(f, t->lz);
		else {
			t->lz = f;
			t->has_lz = true;
		}
	}

	static void apply_key_shift(Node* t, const Key& d) {
		if (!t) return;
		t->key += d;
		t->val = shift_mapping(d, t->val);
		t->prod = shift_mapping(d, t->prod);
		if (t->has_key_lz) t->key_lz += d;
		else {
			t->key_lz = d;
			t->has_key_lz = true;
		}
	}

	static void push(Node* t) {
		if (!t) return;

		if (t->has_lz) {
			apply_val(t->l, t->lz);
			apply_val(t->r, t->lz);
			t->lz = id();
			t->has_lz = false;
		}
		if (t->has_key_lz) {
			apply_key_shift(t->l, t->key_lz);
			apply_key_shift(t->r, t->key_lz);
			t->key_lz = Key();
			t->has_key_lz = false;
		}
	}

	// a の全 key < b の全 key を仮定
	static Node* merge(Node* a, Node* b) {
		if (!a) return b;
		if (!b) return a;
		push(a);
		push(b);
		if (a->pri > b->pri) {
			a->r = merge(a->r, b);
			pull(a);
			return a;
		}
		else {
			b->l = merge(a, b->l);
			pull(b);
			return b;
		}
	}

	// (< key, >= key)
	static pair<Node*, Node*> split_lt(Node* t, const Key& key) {
		if (!t) return { nullptr, nullptr };
		push(t);
		if (t->key < key) {
			auto [a, b] = split_lt(t->r, key);
			t->r = a;
			pull(t);
			return { t, b };
		}
		else {
			auto [a, b] = split_lt(t->l, key);
			t->l = b;
			pull(t);
			return { a, t };
		}
	}

	// (<= key, > key)
	static pair<Node*, Node*> split_le(Node* t, const Key& key) {
		if (!t) return { nullptr, nullptr };
		push(t);
		if (key < t->key) {
			auto [a, b] = split_le(t->l, key);
			t->l = b;
			pull(t);
			return { a, t };
		}
		else {
			auto [a, b] = split_le(t->r, key);
			t->r = a;
			pull(t);
			return { t, b };
		}
	}

	// (< key, = key, > key)
	static tuple<Node*, Node*, Node*> split_eq(Node* t, const Key& key) {
		auto [a, bc] = split_lt(t, key);
		auto [b, c] = split_le(bc, key);
		return { a, b, c };
	}

	static Node* find_node(Node* t, const Key& key) {
		while (t) {
			push(t);
			if (key < t->key) t = t->l;
			else if (t->key < key) t = t->r;
			else return t;
		}
		return nullptr;
	}

	static void clear(Node* t) {
		if (!t) return;
		push(t);
		clear(t->l);
		clear(t->r);
		delete t;
	}

	template<class G>
	static optional<Key> max_right_fail_dfs(Node* t, S& sm, G g) {
		if (!t) return nullopt;
		push(t);

		S nxt = op(sm, prod_all(t->l));
		if (!g(nxt)) return max_right_fail_dfs(t->l, sm, g);
		sm = nxt;

		nxt = op(sm, t->val);
		if (!g(nxt)) return t->key;
		sm = nxt;

		return max_right_fail_dfs(t->r, sm, g);
	}

	template<class G>
	static optional<Key> min_left_fail_dfs(Node* t, S& sm, G g) {
		if (!t) return nullopt;
		push(t);

		S nxt = op(prod_all(t->r), sm);
		if (!g(nxt)) return min_left_fail_dfs(t->r, sm, g);
		sm = nxt;

		nxt = op(t->val, sm);
		if (!g(nxt)) return t->key;
		sm = nxt;

		return min_left_fail_dfs(t->l, sm, g);
	}

	lazy_shift_index_treap() = default;

	lazy_shift_index_treap(const lazy_shift_index_treap&) = delete;
	lazy_shift_index_treap& operator=(const lazy_shift_index_treap&) = delete;

	~lazy_shift_index_treap() {
		clear(root);
	}

	int size() const {
		return node_size(root);
	}

	bool empty() const {
		return root == nullptr;
	}

	bool contains(const Key& key) {
		return find_node(root, key) != nullptr;
	}

	S get(const Key& key) {
		Node* t = find_node(root, key);
		return t ? t->val : e();
	}

	void set(const Key& key, const S& x) {
		auto [a, b, c] = split_eq(root, key);
		clear(b);
		root = merge(a, merge(new Node(key, x, rng()), c));
	}

	void erase(const Key& key) {
		auto [a, b, c] = split_eq(root, key);
		clear(b);
		root = merge(a, c);
	}

	S prod(const Key& l, const Key& r) {
		auto [a, bc] = split_lt(root, l);
		auto [b, c] = split_lt(bc, r);
		S res = prod_all(b);
		root = merge(a, merge(b, c));
		return res;
	}

	S all_prod() const {
		return prod_all(root);
	}

	void apply(const Key& key, const F& f) {
		auto [a, b, c] = split_eq(root, key);
		apply_val(b, f);
		root = merge(a, merge(b, c));
	}

	void apply(const Key& l, const Key& r, const F& f) {
		auto [a, bc] = split_lt(root, l);
		auto [b, c] = split_lt(bc, r);
		apply_val(b, f);
		root = merge(a, merge(b, c));
	}

	void add_ge_index(const Key& x, const Key& d = Key(1)) {
		auto [a, b] = split_lt(root, x);
		apply_key_shift(b, d);
		root = merge(a, b);
	}

	template<class G>
	optional<Key> max_right(const Key& l, G g) {
		assert(g(e()));
		auto [a, b] = split_lt(root, l);
		S sm = e();
		auto res = max_right_fail_dfs(b, sm, g);
		root = merge(a, b);
		return res;
	}

	template<class G>
	optional<Key> min_left(const Key& r, G g) {
		assert(g(e()));
		auto [a, b] = split_lt(root, r);
		S sm = e();
		auto res = min_left_fail_dfs(a, sm, g);
		root = merge(a, b);
		return res;
	}
};

// implicit treap

template<class S, auto op, auto e, class F, auto mapping, auto composition, auto id, auto ts> struct lazy_reversible_treap {
private:
	struct Node {
		S val, prod;
		F lz;
		bool rev;
		uint64_t pri;
		int sz;
		Node* l, * r;
		Node(const S& x, uint64_t p)
			: val(x), prod(x), lz(id()), rev(false), pri(p), sz(1), l(nullptr), r(nullptr) {}
	};

	Node* root = nullptr;

	static uint64_t rng64() {
		static std::mt19937_64 rng(
			(uint64_t)chrono::steady_clock::now().time_since_epoch().count()
		);
		return rng();
	}

	static int size(Node* t) {
		return t ? t->sz : 0;
	}

	static S prod_all(Node* t) {
		return t ? t->prod : e();
	}

	static void all_apply(Node* t, const F& f) {
		if (!t) return;
		t->val = mapping(f, t->val);
		t->prod = mapping(f, t->prod);
		t->lz = composition(f, t->lz);
	}

	static void all_toggle(Node* t) {
		if (!t) return;
		swap(t->l, t->r);
		t->val = ts(t->val);
		t->prod = ts(t->prod);
		t->rev ^= 1;
	}

	static void pull(Node* t) {
		if (!t) return;
		t->sz = 1 + size(t->l) + size(t->r);
		t->prod = op(prod_all(t->l), op(t->val, prod_all(t->r)));
	}

	static void push(Node* t) {
		if (!t) return;
		if (t->rev) {
			all_toggle(t->l);
			all_toggle(t->r);
			t->rev = false;
		}
		if (!(t->lz == id())) {
			all_apply(t->l, t->lz);
			all_apply(t->r, t->lz);
			t->lz = id();
		}
	}

	static Node* merge(Node* a, Node* b) {
		if (!a) return b;
		if (!b) return a;
		if (a->pri > b->pri) {
			push(a);
			a->r = merge(a->r, b);
			pull(a);
			return a;
		}
		else {
			push(b);
			b->l = merge(a, b->l);
			pull(b);
			return b;
		}
	}

	// 左に k 個
	static pair<Node*, Node*> split(Node* t, int k) {
		if (!t) return { nullptr, nullptr };
		push(t);
		if (k <= size(t->l)) {
			auto [a, b] = split(t->l, k);
			t->l = b;
			pull(t);
			return { a, t };
		}
		else {
			auto [a, b] = split(t->r, k - size(t->l) - 1);
			t->r = a;
			pull(t);
			return { t, b };
		}
	}

	static Node* build_rec(const vector<S>& v, int l, int r) {
		Node* t = nullptr;
		for (int i = l; i < r; i++) {
			t = merge(t, new Node(v[i], rng64()));
		}
		return t;
	}

	static void clear(Node* t) {
		if (!t) return;
		push(t);
		clear(t->l);
		clear(t->r);
		delete t;
	}

	static S get_kth(Node* t, int k) {
		while (true) {
			push(t);
			int ls = size(t->l);
			if (k < ls) {
				t = t->l;
			}
			else if (k == ls) {
				return t->val;
			}
			else {
				k -= ls + 1;
				t = t->r;
			}
		}
	}

	static void set_kth(Node* t, int k, const S& x) {
		while (true) {
			push(t);
			int ls = size(t->l);
			if (k < ls) {
				t = t->l;
			}
			else if (k == ls) {
				t->val = x;
				pull(t);
				return;
			}
			else {
				k -= ls + 1;
				t = t->r;
			}
		}
	}

	template <class G>
	static int max_right_dfs(Node* t, S& sm, G g) {
		if (!t) return 0;
		push(t);

		S nxt = op(sm, prod_all(t->l));
		if (!g(nxt)) {
			return max_right_dfs(t->l, sm, g);
		}
		sm = nxt;

		nxt = op(sm, t->val);
		if (!g(nxt)) {
			return size(t->l);
		}
		sm = nxt;

		int res = max_right_dfs(t->r, sm, g);
		return size(t->l) + 1 + res;
	}

	template <class G>
	static int min_left_dfs(Node* t, S& sm, G g) {
		if (!t) return 0;
		push(t);

		S nxt = op(prod_all(t->r), sm);
		if (!g(nxt)) {
			return min_left_dfs(t->r, sm, g);
		}
		sm = nxt;

		nxt = op(t->val, sm);
		if (!g(nxt)) {
			return size(t->r);
		}
		sm = nxt;

		int res = min_left_dfs(t->l, sm, g);
		return size(t->r) + 1 + res;
	}

public:
	lazy_reversible_treap() = default;

	explicit lazy_reversible_treap(int n) {
		vector<S> v(n, e());
		root = build_rec(v, 0, n);
	}

	explicit lazy_reversible_treap(const vector<S>& v) {
		root = build_rec(v, 0, (int)v.size());
	}

	~lazy_reversible_treap() {
		clear(root);
	}

	int size() const {
		return lazy_reversible_treap::size(root);
	}

	bool empty() const {
		return root == nullptr;
	}

	void insert(int p, const S& x) {
		assert(0 <= p && p <= size());
		auto [a, b] = split(root, p);
		root = merge(a, merge(new Node(x, rng64()), b));
	}

	void erase(int p) {
		assert(0 <= p && p < size());
		auto [a, b] = split(root, p);
		auto [c, d] = split(b, 1);
		clear(c);
		root = merge(a, d);
	}

	void set(int p, const S& x) {
		assert(0 <= p && p < size());
		auto [a, b] = split(root, p);
		auto [c, d] = split(b, 1);
		c->val = x;
		pull(c);
		root = merge(a, merge(c, d));
	}

	S get(int p) {
		assert(0 <= p && p < size());
		auto [a, b] = split(root, p);
		auto [c, d] = split(b, 1);
		S res = c->val;
		root = merge(a, merge(c, d));
		return res;
	}

	S prod(int l, int r) {
		assert(0 <= l && l <= r && r <= size());
		auto [a, b] = split(root, l);
		auto [c, d] = split(b, r - l);
		S res = prod_all(c);
		root = merge(a, merge(c, d));
		return res;
	}

	S all_prod() const {
		return prod_all(root);
	}

	void apply(int p, const F& f) {
		assert(0 <= p && p < size());
		auto [a, b] = split(root, p);
		auto [c, d] = split(b, 1);
		all_apply(c, f);
		root = merge(a, merge(c, d));
	}

	void apply(int l, int r, const F& f) {
		assert(0 <= l && l <= r && r <= size());
		auto [a, b] = split(root, l);
		auto [c, d] = split(b, r - l);
		all_apply(c, f);
		root = merge(a, merge(c, d));
	}

	void reverse(int l, int r) {
		assert(0 <= l && l <= r && r <= size());
		auto [a, b] = split(root, l);
		auto [c, d] = split(b, r - l);
		all_toggle(c);
		root = merge(a, merge(c, d));
	}

	template <class G> int max_right(int l, G g) {
		assert(0 <= l && l <= size());
		assert(g(e()));
		auto [a, b] = split(root, l);
		S sm = e();
		int add = max_right_dfs(b, sm, g);
		root = merge(a, b);
		return l + add;
	}

	template <class G> int min_left(int r, G g) {
		assert(0 <= r && r <= size());
		assert(g(e()));
		auto [a, b] = split(root, r);
		S sm = e();
		int add = min_left_dfs(a, sm, g);
		root = merge(a, b);
		return r - add;
	}
};