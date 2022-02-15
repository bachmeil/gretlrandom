import gretl.random;
import std.stdio;

void main() {
	randInit();
	
	foreach(_; 0..20) {
		writeln(shuffleIndex([1, 2, 3, 4, 5]));
	}
	
	randFree();
}
