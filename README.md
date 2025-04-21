# GAMMA_cpp

c++ version of GAMMA(made by Jong Wha J. Joo)<br> 
Original GAMMA paper -> https://doi.org/10.1534/genetics.116.189712

### Prerequisites

1. GAMMA_cpp
```
git clone https://github.com/DGU-CBLAB/Gamma_cpp.git
```

2. Eigen library
```
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
```
```
tar -zxvf ./eigen-3.4.0.tar.gz
```

3. boost library
```
wget https://archives.boost.io/release/1.86.0/source/boost_1_86_0.tar.gz
```
```
tar -zxvf ./boost_1_86_0.tar.gz
```

- Compile
```
cd ./Gamma_cpp
```
```
g++ -O2 -DNDEBUG -pthread -std=c++14 -I ../eigen-3.4.0/ -I ../boost_1_86_0/ ./Gamma_main.cpp ./CBLAB_method.cpp -o Gamma_cpp
```
4. Run
   Gamm_cpp -x x_file_path -y y_file_path -o output_path -t thread_number -p permutation_number 
```bash 
./Gamma_cpp -x ./sample_data/X.txt -y ./sample_data/Y.txt -o ./result.txt -t 1 -p 4
```

5. Result
```
vim ./result.txt
```
