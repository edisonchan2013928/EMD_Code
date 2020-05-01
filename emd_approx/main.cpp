#include "emd_approx.h"

int main(int argc,char**argv)
{
	emd_struct my_emd;

	my_emd.dim = atoi(argv[1]);
	my_emd.method = atoi(argv[2]);
	char*dataFileName = argv[3];
	char*groundFileName = argv[4];
	char*pairFileName = argv[5];
	char*resultFileName = argv[6];
	
	if (my_emd.method == 3)
		my_emd.dim_Red = atoi(argv[7]);
	if (my_emd.method == 4)
		my_emd.hilbert_FileName = argv[7];
	
	//loadData(dataFileName, pairFileName, true, my_emd);
	run_emd_approx(dataFileName, pairFileName, groundFileName, resultFileName, my_emd);
}
