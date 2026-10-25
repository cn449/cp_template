auto centroid_decomposition = [&]() {

	vector<int> subtree_size(n);
	vector<bool> used(n);

	// recalculate size
	auto dfs1 = [&](auto dfs1, int v, int p) -> void {
		subtree_size[v] = 1;
		for (int i : G[v]) {
			if (i == p or used[i]) continue;
			dfs1(dfs1, i, v);
			subtree_size[v] += subtree_size[i];
		}
	};

	// find center
	// cur_sz : 現時点での size
	auto dfs2 = [&](auto dfs2, int root, int pare, int cur_sz) -> void {
		int nxt = -1;
		for (int i : G[root]) {
			if (i == pare or used[i]) continue;
			if (2 * subtree_size[i] > cur_sz) nxt = i;
		}

		if (nxt == -1) {
			// root : 重心
			used[root] = true;

			// ここは ad-hoc に色々やる
			// used なところにはもう触れない

			

			// ここまで

			for (int i : G[root]) {
				if (!used[i]) {
					dfs1(dfs1, i, root);
					dfs2(dfs2, i, root, subtree_size[i]);
				}
			}
		}
		else dfs2(dfs2, nxt, root, cur_sz);
	};

	dfs1(dfs1, 0, -1);
	dfs2(dfs2, 0, -1, n);
};