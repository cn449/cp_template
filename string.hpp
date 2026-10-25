// abababaca -> 123432121
// i を中心とする回文の長さは 2*res[i]-1
vector<int> manacher(const string& s) {
	vector<int> res(ssize(s));
	int i = 0, j = 0;
	while (i < ssize(s)) {
		while (i - j >= 0 && i + j < ssize(s) && s[i - j] == s[i + j]) j++;
		res[i] = j;
		int k = 1;
		while (i - k >= 0 && k + res[i - k] < j) res[i + k] = res[i - k], k++;
		i += k; j -= k;
	}
	return res;
}

template<class T> vector<int> manacher(const vector<T>& s) {
	vector<int> res(ssize(s));
	int i = 0, j = 0;
	while (i < ssize(s)) {
		while (i - j >= 0 && i + j < ssize(s) && s[i - j] == s[i + j]) j++;
		res[i] = j;
		int k = 1;
		while (i - k >= 0 && k + res[i - k] < j) res[i + k] = res[i - k], k++;
		i += k; j -= k;
	}
	return res;
}

// s, t に $ が含まれないことを仮定している
// (s の開始 index, t の開始 index, 長さ) を返す
array<int, 3> longest_common_substring(const string& s, const string& t) {
	string w = s + "$" + t;
	auto sa = suffix_array(w);
	auto lcp = lcp_array(w, sa);
	int len = -1, idx = 0;
	for (int i = 1; i < ssize(s) + ssize(t); i++) {
		bool f0 = sa[i] < ssize(s);
		bool f1 = sa[i + 1] < ssize(s);
		if (f0 != f1) if (chmax(len, lcp[i])) idx = i;
	}
	array<int, 3> res = { sa[idx], sa[idx + 1], len };
	if (res[0] > res[1]) swap(res[0], res[1]);
	res[1] -= ssize(s) + 1;
	return res;
}

vector<int> Z_algorithm(string s) {
	int n = s.size();
	vector<int> z(n);
	z[0] = n;
	int i = 1, j = 0;
	while (i < n) {
		while (i + j < n && s[j] == s[i + j]) j++;
		z[i] = j;
		if (j == 0) {
			i++; continue;
		}
		int k = 1;
		while (k < j && k + z[k] < j) {
			z[i + k] = z[k]; k++;
		}
		i += k; j -= k;
	}
	return z;
}

// sort substring of s
vector<P<int>> sort_substring(string s) {
	int n = s.size();
	vector<vector<P<int>>> v(26);
	rep(i, n) v[s[i] - 'a'].push_back({ i, i + 1 });
	vector<P<int>> res;

	auto dfs = [&](auto dfs, vector<P<int>> a) -> void {
		if (a.empty()) return;
		for (auto p : a) res.push_back(p);
		vector<vector<P<int>>> w(26);
		for (auto p : a) {
			if (p.second != n) {
				char c = s[p.second];
				w[c - 'a'].push_back({ p.first, p.second + 1 });
			}
		}
		rep(i, 26) if (!w[i].empty()) dfs(dfs, w[i]);
	};

	rep(i, 26) dfs(dfs, v[i]);
	return res;
}

// sort substring of s, unique
vector<P<int>> sort_substring_unique(string s) {
	int n = s.size();
	vector<vector<P<int>>> v(26);
	rep(i, n) v[s[i] - 'a'].push_back({ i, i + 1 });
	vector<P<int>> res;

	auto dfs = [&](auto dfs, vector<P<int>> a) -> void {
		if (a.empty()) return;
		res.push_back(a[0]);
		vector<vector<P<int>>> w(26);
		for (auto p : a) {
			if (p.second != n) {
				char c = s[p.second];
				w[c - 'a'].push_back({ p.first, p.second + 1 });
			}
		}
		rep(i, 26) if (!w[i].empty()) dfs(dfs, w[i]);
	};

	rep(i, 26) dfs(dfs, v[i]);
	return res;
}

constexpr int SIZE = 26;

class AhoCorasick {
	struct Node {
		int par = -1;
		int go[SIZE];
		int failure_link = 0;
		bool last = false;
		Node() { fill(go, go + SIZE, -1); }
	};
public:
	vector<Node> trie;

	AhoCorasick(vector<string>& s) {
		trie.resize(1);
		for (auto& t : s) {
			int cur = 0;
			for (char c : t) {
				int v = c - 'a';
				if (trie[cur].go[v] == -1) {
					trie[cur].go[v] = ssize(trie);
					trie.emplace_back();
					trie.back().par = cur;
				}
				cur = trie[cur].go[v];
			}
			trie[cur].last = true;
		}
		queue<int> q;
		q.push(0);
		while (!q.empty()) {
			int v = q.front(); q.pop();
			trie[v].last |= trie[trie[v].failure_link].last;
			rep(i, SIZE) {
				int w = trie[v].go[i];
				if (w == -1) continue;
				if (v != 0) {
					int f = trie[v].failure_link;
					while (f && trie[f].go[i] == -1) f = trie[f].failure_link;
					if (trie[f].go[i] != -1) trie[w].failure_link = trie[f].go[i];
				}
				q.push(w);
			}
		}
	}

	int size() { return ssize(trie); }
	int parent(int v) { return trie[v].par; }
	int next_inner(int v, char c) { return trie[v].go[c - 'a']; }
	int failure(int v) { return trie[v].failure_link; }
	int next(int v, char c) {
		while (v != 0 && next_inner(v, c) == -1) v = failure(v);
		v = next_inner(v, c);
		if (v == -1) v = 0;
		return v;
	}
	bool is_last(int v) { return trie[v].last; }

};