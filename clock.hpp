class Timer {
private:
	chrono::high_resolution_clock::time_point start_time;
public:
	void start() { start_time = chrono::high_resolution_clock::now(); }
	ll elapsed_ms() const {
		auto now = chrono::high_resolution_clock::now();
		return chrono::duration_cast<chrono::milliseconds>(now - start_time).count();
	}
};