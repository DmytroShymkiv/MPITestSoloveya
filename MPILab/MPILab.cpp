#include <iostream>
#include <math.h>
#include <mpi.h>
#include <cstdlib>
#include "BigInt.h"

using namespace std;

long long binpow(long long a, long long e, long long m) {
	long long res = 1;
	while (e != 0) {
		if (e % 2 == 1)
			res = (res * a) % m;

		e = e / 2;
		a = (a * a) % m;

	}
	return res % m;
}
long long gcd(long long a, long long m)
{
	while (a != 0 && m != 0) {
		if (a > m)
			a = a % m;
		else
			m = m % a;
	}
	return a + m;
}
int JacobiSymbol(long long a1, long long p1)
{
	if (p1 < 1 || (p1 % 2 == 0))
		return 10;
	if (gcd(a1, p1) != 1)
		return 0;
	long long a = a1;
	long long p = p1;
	a = a % p;
	long long t = 1;
	while (a != 0) {
		while ((a % 2 == 0)) {
			a = a / 2;
			long long r = p % 8;
			if (r == 3 || r == 5)
				t = -t;
		}
		long long x = p;
		p = a;
		a = x;
		if (a % 4 == 3 && p % 4 == 3)
			t = -t;

		a = a % p;
	}
	if (p == 1)
		return t;
	return 0;
}
int TestSoloveya(long long a,int N)
{
	int res = 1;
	int i;

	for (int i = 0; i < N; i += 1)
	{
		long long k = 1 + rand() % (a - 1);

		if (gcd(a, k) > 1)
		{
			res = 0;
			break;
		}

		long long j = binpow(k, ((a - 1) / 2), a);
		int jacobiSymbol = JacobiSymbol(k, a);
		if (jacobiSymbol < 0)
			j = j - a;

		if (j != jacobiSymbol)
		{
			res = 0;
			break;
		}

	}
	return res;
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank;
    int size;

	double start, end;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	long number;
	int N;
	if (rank == 0) {
		N = 10000000 / size;
		number = atoi(argv[1]);
	}
		

	MPI_Bcast(&number, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	start = MPI_Wtime();

	int result = TestSoloveya(number,N); 
	int globalRes;
	MPI_Reduce(&result, &globalRes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&globalRes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	end = MPI_Wtime();

	if (globalRes > 0)
		globalRes = 1;

	if (rank == 0) {
		cout << "Result " << globalRes << endl;
		cout << "Time " << (end - start) << endl;
	}
    MPI_Finalize();
	return 0;
}
