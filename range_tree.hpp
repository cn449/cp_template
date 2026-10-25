// https://kopricky.github.io/code/SegmentTrees/rangetree_pointupdate.html
// ↑ をもとに自分が使いやすい形に改良

template<typename CandidateType, typename ValueType, auto op, auto e> class RangeTree {
public:
	static_assert(is_integral<CandidateType>::value, "Integral required.");
private:
	using CT = CandidateType;
	using VT = ValueType;
	using pcc = pair<CT, CT>;
	using pci = pair<CT, int>;
	int n, sz;
	vector<segtree<VT, op, e>> seg;
	vector<vector<pcc>> yx;
	vector<pcc> sorted;
	void set_(int id, const CT x, const CT y, const VT val) {
		id += n - 1;
		const int yid = lower_bound(all(yx[id]), pcc(y, x)) - yx[id].begin();
		seg[id].set(yid, val);
		while (id > 0) {
			id = (id - 1) >> 1;
			const int yid = lower_bound(all(yx[id]), pcc(y, x)) - yx[id].begin();
			seg[id].set(yid, val);
		}
	}
	VT prod(const int lxid, const int rxid, const CT ly, const CT ry, const int k, const int l, const int r) {
		if (r <= lxid || rxid <= l) return e();
		if (lxid <= l && r <= rxid) {
			const int lyid = lower_bound(all(yx[k]), pcc(ly, numeric_limits<CT>::min())) - yx[k].begin();
			const int ryid = upper_bound(all(yx[k]), pcc(ry, numeric_limits<CT>::min())) - yx[k].begin();
			return (lyid >= ryid) ? e() : seg[k].prod(lyid, ryid);
		}
		else {
			int m = midpoint(l, r);
			return op(prod(lxid, rxid, ly, ry, 2 * k + 1, l, m), prod(lxid, rxid, ly, ry, 2 * k + 2, m, r));
		}
	}
public:
	RangeTree(const vector<pcc>& cand, const vector<VT>& val) : n(1), sz(ssize(cand)), sorted(sz) {
		n = bit_ceil(uint(sz));
		rep(i, sz) sorted[i] = { cand[i].first, i };
		sort(all(sorted), [&](const pcc& a, const pcc& b) { return (a.first == b.first) ? (cand[a.second].second < cand[b.second].second) : (a.first < b.first); });
		yx.resize(2 * n - 1), seg.resize(2 * n - 1);
		rep(i, sz) {
			yx[i + n - 1] = { {sorted[i].second, sorted[i].first} };
			vector<VT> arg = { val[sorted[i].second] };
			seg[i + n - 1] = segtree<VT, op, e>(arg);
			sorted[i].second = cand[sorted[i].second].second;
		}
		for (int i = n - 2; i >= 0; --i) {
			yx[i].resize(ssize(yx[2 * i + 1]) + ssize(yx[2 * i + 2]));
			if (yx[i].empty()) continue;
			merge(all(yx[2 * i + 1]), all(yx[2 * i + 2]), yx[i].begin(), [&](const pcc& a, const pcc& b) { return (cand[a.first].second == cand[b.first].second) ? (a.second < b.second) : (cand[a.first].second < cand[b.first].second); });
			vector<VT> arg(ssize(yx[i]));
			rep(j, ssize(yx[i])) arg[j] = val[yx[i][j].first];
			seg[i] = segtree<VT, op, e>(arg);
		}
		rep(i, 2 * n - 1) for (pcc& e0 : yx[i]) e0.first = cand[e0.first].second;
	}
	void set(const CT x, const CT y, const VT val) {
		const int id = lower_bound(all(sorted), pcc(x, y)) - sorted.begin();
		return set_(id, x, y, val);
	}
	VT prod(const CT lx, const CT rx, const CT ly, const CT ry) {
		const int lxid = lower_bound(all(sorted), pcc(lx, numeric_limits<CT>::min())) - sorted.begin();
		const int rxid = upper_bound(all(sorted), pcc(rx, numeric_limits<CT>::min())) - sorted.begin();
		return (lxid >= rxid) ? e() : prod(lxid, rxid, ly, ry, 0, 0, n);
	}
};