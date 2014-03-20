
size_t _Find_first(bitset<MAX_EDGES>& vSet){
	for (int i = 0; i < vSet.size(); i++){
		if (vSet[i])
			return i;
	}
	return vSet.size();
}

size_t _Find_next(bitset<MAX_EDGES>& vSet, size_t v){
	for (int i = v+1; i < vSet.size(); i++){
		if (vSet[i])
			return i;
	}
	return vSet.size();
}

struct timeval{
	int tv_sec;
	int tv_usec;
};

void gettimeofday(timeval* time, int dummy){
	time->tv_sec = 100;
	time->tv_usec = 200;
}
