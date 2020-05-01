g++ -c emd.cpp -w -std=c++11
g++ -c emd_approx.cpp -w -std=c++11
g++ -c exact_emd.cpp -w -std=c++11
g++ -c init.cpp -w -std=c++11
g++ -c LB_dimRed.cpp -w -std=c++11
g++ -c LB_IM.cpp -w -std=c++11
g++ -c LB_Proj.cpp -w -std=c++11
g++ -c UB_G.cpp -w -std=c++11
g++ -c UB_H.cpp -w -std=c++11

g++ main.cpp -w -O3 -o main emd.o emd_approx.o exact_emd.o init.o LB_dimRed.o LB_IM.o LB_Proj.o UB_G.o UB_H.o

#my_emd.dim = atoi(argv[1]);
#my_emd.method = atoi(argv[2]);
#char*dataFileName = argv[3];
#char*groundFileName = argv[4];
#char*pairFileName = argv[5];
#char*resultFileName = argv[6];
#if (my_emd.method == 3)
#	my_emd.dim_Red = atoi(argv[7]);
#if (my_emd.method == 4)
#	my_emd.hilbert_FileName = argv[7];

#method = 0 (exact EMD)
#method = 1 (LB_Proj)
#method = 2 (LB_IM)
#method = 3 (LB_Red)
#method = 4 (UB_H)
#method = 5 (UB_G)

dim=64
method=1 #LB_Proj
dataFileName="./Datasets/256_ObjectCategories_64RGB_FV/256_ObjectCategories_64RGB_FVMatrix.txt"
groundFileName="./Datasets/ground_vector_64RGB.txt"
pairFileName="./Datasets/Pair_File/256_OC_1000Pairs.txt"
resultFileName="./Result/256_OC_M1.txt"

./main $dim $method $dataFileName $groundFileName $pairFileName $resultFileName

#Output:
#Method = 1: (throughput (queries/sec))
