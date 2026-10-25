template<class S, auto op, auto e> class Dynamic_Segtree {
	struct Node { S val = e(); Node* l = nullptr, * r = nullptr, * par = nullptr; };
	ll n; int d;
	ll _n;
	Node* root;
	S val_of_node(Node* s) {
		if (s == nullptr) return e();
		else return s->val;
	}

	template<class G> ll max_right(Node* node, ll a, ll b, ll l, ll r, G& g, S& sm) {
		if (b <= l || r <= a) return r;
		if (l <= a && b <= r) {
			S nxt = op(sm, val_of_node(node));
			if (g(nxt)) {
				sm = nxt;
				return r;
			}
			if (b - a == 1) return a;
		}
		ll m = midpoint(a, b);
		ll res = max_right(node == nullptr ? nullptr : node->l, a, m, l, r, g, sm);
		if (res != r) return res;
		return max_right(node == nullptr ? nullptr : node->r, m, b, l, r, g, sm);
	}
	template<class G> ll min_left(Node* node, ll a, ll b, ll l, ll r, G& g, S& sm) {
		if (b <= l || r <= a) return l;
		if (l <= a && b <= r) {
			S nxt = op(val_of_node(node), sm);
			if (g(nxt)) {
				sm = nxt;
				return l;
			}
			if (b - a == 1) return b;
		}
		ll m = midpoint(a, b);
		ll res = min_left(node == nullptr ? nullptr : node->r, m, b, l, r, g, sm);
		if (res != l) return res;
		return min_left(node == nullptr ? nullptr : node->l, a, m, l, r, g, sm);
	}
public:
	Dynamic_Segtree(ll n_) {
		_n = n_;
		n = 1, d = 0;
		while (n < n_) { n <<= 1; d++; }
		root = new Node();
	}
	void set(ll pos, S val) {
		Node* cur = root;
		for (int i = d - 1; i >= 0; i--) {
			if (!(pos >> i & 1)) {
				if (!cur->l) {
					cur->l = new Node();
					cur->l->par = cur;
				}
				cur = cur->l;
			}
			else {
				if (!cur->r) {
					cur->r = new Node();
					cur->r->par = cur;
				}
				cur = cur->r;
			}
		}
		cur->val = val;
		while (cur->par != nullptr) {
			cur = cur->par;
			cur->val = op(val_of_node(cur->l), val_of_node(cur->r));
		}
	}
	S prod(Node* node, ll a, ll b, ll l, ll r) {
		if (node == nullptr or r <= a or b <= l) return e();
		if (l <= a && b <= r) return node->val;
		ll m = midpoint(a, b);
		return op(prod(node->l, a, m, l, r), prod(node->r, m, b, l, r));
	}
	S prod(ll l, ll r) { return prod(root, 0, n, l, r); }
	S all_prod() { return val_of_node(root); }
	S get(ll p) { return prod(p, p + 1); }

	template<bool (*g)(S)> ll max_right(ll l) { return max_right(l, [](S x) { return g(x); }); }
	template<class G> ll max_right(ll l, G g) {
		assert(0 <= l && l <= _n);
		assert(g(e()));
		if (l == _n) return _n;
		S sm = e();
		return max_right(root, 0, n, l, _n, g, sm);
	}
	template<bool (*g)(S)> ll min_left(ll r) { return min_left(r, [](S x) { return g(x); }); }
	template<class G> ll min_left(ll r, G g) {
		assert(0 <= r && r <= _n);
		assert(g(e()));
		if (r == 0) return 0;
		S sm = e();
		return min_left(root, 0, n, 0, r, g, sm);
	}
};

template<class S, auto op, auto e, class F, auto mapping, auto composition, auto id, auto init> class Lazy_Dynamic_Segtree {
	struct Node {
		S val;
		F laz;
		bool has_lazy;
		Node* l;
		Node* r;
		Node* par;
		Node(const S& v, Node* p = nullptr) : val(v), laz(id()), has_lazy(false), l(nullptr), r(nullptr), par(p) {}
	};
	ll n; int d;
	ll _n;
	Node* root;
	S get_val(Node* node, ll l, ll r) const {
		if (node == nullptr) return init(l, r);
		return node->val;
	}
	Node* make_node(ll l, ll r, Node* par = nullptr) { return new Node(init(l, r), par); }
	void all_apply(Node* node, const F& f) {
		node->val = mapping(f, node->val);
		if (node->has_lazy) node->laz = composition(f, node->laz);
		else {
			node->laz = f;
			node->has_lazy = true;
		}
	}
	void push(Node* node, ll l, ll r) {
		if (node == nullptr) return;
		if (!node->has_lazy) return;
		if (r - l == 1) {
			node->laz = id();
			node->has_lazy = false;
			return;
		}
		ll m = midpoint(l, r);
		if (!node->l) node->l = make_node(l, m, node);
		if (!node->r) node->r = make_node(m, r, node);
		all_apply(node->l, node->laz);
		all_apply(node->r, node->laz);
		node->laz = id();
		node->has_lazy = false;
	}
	void update(Node* node, ll l, ll r) {
		ll m = midpoint(l, r);
		node->val = op(get_val(node->l, l, m), get_val(node->r, m, r));
	}
	void set(Node* node, ll l, ll r, ll p, const S& x) {
		if (r - l == 1) {
			node->val = x;
			node->laz = id();
			node->has_lazy = false;
			return;
		}
		push(node, l, r);
		ll m = midpoint(l, r);
		if (p < m) {
			if (!node->l) node->l = make_node(l, m, node);
			set(node->l, l, m, p, x);
		}
		else {
			if (!node->r) node->r = make_node(m, r, node);
			set(node->r, m, r, p, x);
		}
		update(node, l, r);
	}
	S prod(Node* node, ll a, ll b, ll l, ll r) {
		if (r <= a || b <= l) return e();
		if (node == nullptr) {
			if (l <= a && b <= r) return init(a, b);
			ll m = midpoint(a, b);
			return op(prod(nullptr, a, m, l, r), prod(nullptr, m, b, l, r));
		}
		if (l <= a && b <= r) return node->val;
		push(node, a, b);
		ll m = midpoint(a, b);
		return op(prod(node->l, a, m, l, r), prod(node->r, m, b, l, r));
	}
	void apply(Node*& node, Node* par, ll a, ll b, ll l, ll r, const F& f) {
		if (r <= a || b <= l) return;
		if (node == nullptr) node = make_node(a, b, par);
		if (l <= a && b <= r) { all_apply(node, f); return; }
		push(node, a, b);
		ll m = midpoint(a, b);
		apply(node->l, node, a, m, l, r, f);
		apply(node->r, node, m, b, l, r, f);
		update(node, a, b);
	}

	template<class G> ll max_right(Node* node, ll a, ll b, ll l, ll r, G& g, S& sm) {
		if (b <= l || r <= a) return r;
		if (l <= a && b <= r) {
			S nxt = op(sm, get_val(node, a, b));
			if (g(nxt)) {
				sm = nxt;
				return r;
			}
			if (b - a == 1) return a;
		}
		if (node != nullptr) push(node, a, b);
		ll m = midpoint(a, b);
		ll res = max_right(node == nullptr ? nullptr : node->l, a, m, l, r, g, sm);
		if (res != r) return res;
		return max_right(node == nullptr ? nullptr : node->r, m, b, l, r, g, sm);
	}
	template<class G> ll min_left(Node* node, ll a, ll b, ll l, ll r, G& g, S& sm) {
		if (b <= l || r <= a) return l;
		if (l <= a && b <= r) {
			S nxt = op(get_val(node, a, b), sm);
			if (g(nxt)) {
				sm = nxt;
				return l;
			}
			if (b - a == 1) return b;
		}
		if (node != nullptr) push(node, a, b);
		ll m = midpoint(a, b);
		ll res = min_left(node == nullptr ? nullptr : node->r, m, b, l, r, g, sm);
		if (res != l) return res;
		return min_left(node == nullptr ? nullptr : node->l, a, m, l, r, g, sm);
	}
public:
	Lazy_Dynamic_Segtree(ll n_) {
		_n = n_;
		n = 1, d = 0;
		while (n < n_) { n <<= 1; d++; }
		root = make_node(0, n);
	}
	void set(ll p, const S& x) { set(root, 0, n, p, x); }
	S prod(ll l, ll r) { return prod(root, 0, n, l, r); }
	// S all_prod() const { return root->val; }
	S get(ll p) { return prod(p, p + 1); }
	void apply(ll l, ll r, const F& f) { apply(root, nullptr, 0, n, l, r, f); }
	void apply(ll p, const F& f) { apply(p, p + 1, f); }

	template<bool (*g)(S)> ll max_right(ll l) { return max_right(l, [](S x) { return g(x); }); }
	template<class G> ll max_right(ll l, G g) {
		assert(0 <= l && l <= _n);
		assert(g(e()));
		if (l == _n) return _n;
		S sm = e();
		return max_right(root, 0, n, l, _n, g, sm);
	}
	template<bool (*g)(S)> ll min_left(ll r) { return min_left(r, [](S x) { return g(x); }); }
	template<class G> ll min_left(ll r, G g) {
		assert(0 <= r && r <= _n);
		assert(g(e()));
		if (r == 0) return 0;
		S sm = e();
		return min_left(root, 0, n, 0, r, g, sm);
	}
};