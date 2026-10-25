istream& operator >>(istream& s, __int128_t& v) {
	string str; s >> str;
	bool neg = false;
	v = 0;
	for (char c : str) {
		if (c == '-') neg = true;
		else { v *= 10; v += c - '0'; }
	}
	if (neg) v *= -1;
	return s;
}

ostream& operator <<(ostream& s, __int128_t v) {
	if (v == 0) { return s << 0; }
	string str;
	bool neg = false;
	if (v < 0) { neg = true; v *= -1; }
	while (v > 0) {
		str += '0' + int(v % 10);
		v /= 10;
	}
	if (neg) str += '-';
	reverse(all(str));
	return s << str;
}
