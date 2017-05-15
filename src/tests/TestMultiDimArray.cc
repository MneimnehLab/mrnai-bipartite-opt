#include <iostream>
using namespace std;

#include "../utils/MultiIndexArray.h"

int main()
{
	// MultiIndexArray<double> mia ({4, 5, 3});

	vector<int> v {4,5,3};

	MultiIndexArray<double> mia;
	cout << "A" << endl;
	mia = MultiIndexArray<double>::CreateMultiIndexArray(v);

	cout << "X" << endl;
	int c = 0;
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<5; j++)	
		{
			for(int k=0; k<3; k++)
			{
				int indices[] = {i,j,k};
				mia[indices] = (double)c++ + 0.5;
				// cout << mia.doHash(indices) << endl;;
				// cout << "B" << endl;
			}
		}
	}

	for(int i=0; i<4*5*3; i++)
	{
		cout << mia[i] << " ";
	}
	cout << endl;

	return 0;
}
