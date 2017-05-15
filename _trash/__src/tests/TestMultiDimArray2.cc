#include <iostream>
using namespace std;

#include "../utils/MultiIndexArray.h"

int main()
{
	// MultiIndexArray<double> mia ({4, 5, 3});

	vector<int> v {8+1,4+1,3+1,5+1};

	MultiIndexArray<double> mia;
	cout << "A" << endl;
	mia = MultiIndexArray<double>::CreateMultiIndexArray(v);

	cout << "X" << endl;
	
	// for(int i=0; i<8*4*3*5; i++)
	for(int i=0; i<9*5*4*6; i++)
	{
		// vector<int> rev = mia.revHash(i);
		int index[4];
		mia.revHash(i, index);
		cout << i << "\t";

		for(int j=0; j<4; j++)
			cout << index[j] << " ";
		
		cout << "\t > ";

		int h = mia.doHash(index);
		cout << h;


		cout << endl;
	}
	cout << endl;

	return 0;
}
