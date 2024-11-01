# GAMMA_cpp

c++ version of GAMMA(made by Jong Wha J. Joo)<br> 
Original GAMMA paper -> https://doi.org/10.1534/genetics.116.189712

### Prerequisites

1. GAMMA_cpp
```
git clone https://github.com/taegun89/Gamma_cpp.git
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
   Gamm_cpp <Genotype file> <Phenotypes file> <threadNums> <output> 
```bash 
./Gamma_cpp ./sample_data/X.txt ./sample_data/Y.txt 1 ./result.txt
```

5. Result
```
vim ./result.txt
```
