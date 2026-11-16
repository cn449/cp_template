struct nimber_info {
	uint8_t table[256][256]{};
	uint8_t inv8[256]{}, norm8[256]{};
	uint16_t mul15[2][256]{}, norm16[2][256]{};
	uint32_t mul31[4][256]{}, norm32[4][256]{};
	ull sq[8][256]{};
	uint16_t prod16(uint16_t x, uint16_t y) const {
		uint8_t a = x, b = x >> 8, c = y, d = y >> 8;
		uint8_t ac = table[a][c], bd = table[b][d];
		return uint16_t(ac ^ table[bd][128]) | (uint16_t(table[a ^ b][c ^ d] ^ ac) << 8);
	}
	uint32_t prod32(uint32_t x, uint32_t y) const {
		uint16_t a = x, b = x >> 16, c = y, d = y >> 16;
		uint16_t ac = prod16(a, c), bd = prod16(b, d);
		uint16_t low = ac ^ mul15[0][bd & 255] ^ mul15[1][bd >> 8];
		return uint32_t(low) | (uint32_t(prod16(a ^ b, c ^ d) ^ ac) << 16);
	}
	ull prod64(ull x, ull y) const {
		uint32_t a = x, b = x >> 32, c = y, d = y >> 32;
		uint32_t ac = prod32(a, c), bd = prod32(b, d);
		uint32_t low = ac ^ mul31[0][bd & 255] ^ mul31[1][(bd >> 8) & 255] ^ mul31[2][(bd >> 16) & 255] ^ mul31[3][bd >> 24];
		return ull(low) | (ull(prod32(a ^ b, c ^ d) ^ ac) << 32);
	}
	uint16_t inverse16(uint16_t x) const {
		if (x < 256) return inv8[x];
		uint8_t a = x, b = x >> 8;
		uint8_t t = inv8[table[a][a ^ b] ^ norm8[b]];
		return uint16_t(table[a ^ b][t]) | (uint16_t(table[b][t]) << 8);
	}
	uint32_t inverse32(uint32_t x) const {
		if (x < 65536) return inverse16(x);
		uint16_t a = x, b = x >> 16;
		uint16_t norm = prod16(a, a ^ b) ^ norm16[0][b & 255] ^ norm16[1][b >> 8];
		uint16_t t = inverse16(norm);
		return uint32_t(prod16(a ^ b, t)) | (uint32_t(prod16(b, t)) << 16);
	}
	ull inverse64(ull x) const {
		if (x <= UINT32_MAX) return inverse32(x);
		uint32_t a = x, b = x >> 32;
		uint32_t norm = prod32(a, a ^ b) ^ norm32[0][b & 255] ^ norm32[1][(b >> 8) & 255] ^ norm32[2][(b >> 16) & 255] ^ norm32[3][b >> 24];
		uint32_t t = inverse32(norm);
		return ull(prod32(a ^ b, t)) | (ull(prod32(b, t)) << 32);
	}
	ull square64(ull x) const {
		return sq[0][x & 255] ^ sq[1][(x >> 8) & 255] ^ sq[2][(x >> 16) & 255] ^ sq[3][(x >> 24) & 255] ^ sq[4][(x >> 32) & 255] ^ sq[5][(x >> 40) & 255] ^ sq[6][(x >> 48) & 255] ^ sq[7][x >> 56];
	}
	nimber_info() {
		table[1][1] = 1;
		for (int h : {1, 2, 4}) {
			int s = 1 << (h * 2), mask = (1 << h) - 1;
			for (int x = 0; x < s; ++x) for (int y = 0; y < s; ++y) {
				int a = x & mask, b = x >> h, c = y & mask, d = y >> h;
				int ac = table[a][c], bd = table[b][d];
				table[x][y] = (ac ^ table[bd][1 << (h - 1)]) | ((table[a ^ b][c ^ d] ^ ac) << h);
			}
		}
		for (int x = 1; x < 256; ++x) for (int y = 1; y < 256; ++y) if (table[x][y] == 1) inv8[x] = y;
		for (int x = 0; x < 256; ++x) norm8[x] = table[table[x][x]][128];
		for (int i = 0; i < 2; ++i) for (int x = 0; x < 256; ++x) {
			uint16_t a = uint16_t(x) << (i * 8);
			mul15[i][x] = prod16(a, uint16_t(1) << 15);
			norm16[i][x] = prod16(prod16(a, a), uint16_t(1) << 15);
		}
		for (int i = 0; i < 4; ++i) for (int x = 0; x < 256; ++x) {
			uint32_t a = uint32_t(x) << (i * 8);
			mul31[i][x] = prod32(a, uint32_t(1) << 31);
			norm32[i][x] = prod32(prod32(a, a), uint32_t(1) << 31);
		}
		for (int i = 0; i < 8; ++i) for (int x = 0; x < 256; ++x) {
			ull a = ull(x) << (i * 8);
			sq[i][x] = prod64(a, a);
		}
	}
};

class nimber {
	inline static const nimber_info info{};
	ull n;
public:
	nimber() : n(0) {}
	template<integral T> nimber(T n_) : n(ull(n_)) {}
	nimber operator -() const { return *this; }
	bool operator == (const nimber& m) const { return n == m.n; }
	bool operator != (const nimber& m) const { return n != m.n; }
	nimber& operator += (const nimber& m) { n ^= m.n; return *this; }
	nimber& operator -= (const nimber& m) { n ^= m.n; return *this; }
	nimber& operator *= (const nimber& m) { n = info.prod64(n, m.n); return *this; }
	nimber& operator /= (const nimber& m) { assert(m.n != 0); n = info.prod64(n, info.inverse64(m.n)); return *this; }
	nimber operator +(nimber m) const { return nimber(*this) += m; }
	nimber operator -(nimber m) const { return nimber(*this) -= m; }
	nimber operator *(nimber m) const { return nimber(*this) *= m; }
	nimber operator /(nimber m) const { return nimber(*this) /= m; }
	nimber& operator ++ () { n ^= 1; return *this; }
	nimber operator ++ (int) { nimber m = *this; n ^= 1; return m; }
	nimber& operator -- () { n ^= 1; return *this; }
	nimber operator -- (int) { nimber m = *this; n ^= 1; return m; }
	nimber inverse() const { assert(n != 0); return info.inverse64(n); }
	nimber square() const { return info.square64(n); }
	nimber pow(ull k) const {
		nimber a = *this, r = 1;
		while (k) {
			if (k & 1) r *= a;
			k >>= 1;
			if (k) a = a.square();
		}
		return r;
	}
	ull val() const { return n; }
	friend istream& operator >> (istream& s, nimber& a) {
		ull x;
		if (s >> x) a = nimber(x);
		return s;
	}
	friend ostream& operator << (ostream& s, const nimber& a) { return s << a.n; }
};