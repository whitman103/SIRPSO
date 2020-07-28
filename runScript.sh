g++ -O2 takeInput.cpp -o takeInput.exe
mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputOne.txt"
mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputTwo.txt"
mpirun -np 20 ./SIRPSO_MPI.exe "20" "InputFolder//outputThree.txt"
