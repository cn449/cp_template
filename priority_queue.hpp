// ABC306-E
template<class T> class del_pql {
	priority_queue<T> que;
	priority_queue<T> rm_que;
public:
	void op() { while (!que.empty() && !rm_que.empty() && que.top() == rm_que.top()) { que.pop(); rm_que.pop(); } }
	void push(T a) { que.push(a); }
	void pop() { op();  que.pop(); }
	T top() { op(); return que.top(); }
	void erase(T a) { rm_que.push(a); }
	bool empty() { return que.size() == rm_que.size(); }
	int size() { return que.size() - rm_que.size(); }
};

template<class T> class del_pqs {
	priority_queue<T, vector<T>, greater<T>> que;
	priority_queue<T, vector<T>, greater<T>> rm_que;
public:
	void op() { while (!que.empty() && !rm_que.empty() && que.top() == rm_que.top()) { que.pop(); rm_que.pop(); } }
	void push(T a) { que.push(a); }
	void pop() { op(); que.pop(); }
	T top() { op(); return que.top(); }
	void erase(T a) { rm_que.push(a); }
	bool empty() { return que.size() == rm_que.size(); }
	int size() { return que.size() - rm_que.size(); }
};

// https://judge.yosupo.jp/submission/269462
template<class T> class Double_Ended_priority_queue {
	del_pql<T> que_large;
	del_pqs<T> que_small;
public:
	void push(T a) { que_large.push(a); que_small.push(a); }
	void pop_large() {
		T a = que_large.top(); que_large.pop();
		que_small.erase(a);
	}
	void pop_small() {
		T a = que_small.top(); que_small.pop();
		que_large.erase(a);
	}
	T get_large() { return que_large.top(); }
	T get_small() { return que_small.top(); }
	bool empty() { return que_large.empty(); }
	int size() { return que_large.size(); }
};