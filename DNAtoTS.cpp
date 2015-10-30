//code for DNA conversion: initial_value=0, A=+2, G=+1, C=-1, T=-2
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

char *INFILENAME; //set the path and name of input DNA file
char *OUTFILENAME; //set the path and name for output Time series;

int main(int argc , char *argv[])
{
	char ncl;
	int val = 0;	//value of TS is initially set to 0, val changes with DNA sequence
	long long cnt = 0;	//count of total number of base pairs in the DNA
	double t1,t2;

	INFILENAME = argv[1];
	OUTFILENAME = argv[2];

	ifstream infile(INFILENAME);	//input DNA file
	ofstream outfile(OUTFILENAME);	//output file for the TS

	t1 = clock();

	if (infile.is_open())
	{
		while(infile.get(ncl) && cnt < 1000000)
		{

			if (ncl=='A' || ncl=='a')
			{
				val = val+2;
				cnt++;
				outfile<<val<<'\n';
			}
			else if (ncl=='G' || ncl=='g')
			{
				val = val+1;
				cnt++;
				outfile<<val<<'\n';
			}
			else if (ncl=='C' || ncl=='c')
			{
				val = val-1;
				cnt++;
				outfile<<val<<'\n';
			}
			else if (ncl=='T' || ncl=='t')
			{
				val = val-2;
				cnt++;
				outfile<<val<<'\n';
			}
		}
	}
	else
		cout<<"Cannot open the DNA file!";




	t2 = clock();
	cout << "Total Execution Time : "<<((t2-t1)/CLOCKS_PER_SEC)<<endl;
	cout << "Total number of nucleotides in the DNA: "<<cnt << endl;

	infile.close();
	outfile.close();

	return 0;
}
