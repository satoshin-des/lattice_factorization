# lattice_factorization

このリポジトリに含まれているソースたちは，大学院での研究の際に作成されたものです．以下，その論文：
<br>[
近似最近ベクトル探索と埋め込み法を用いた格子による素因数分解法の実装報告](https://www.iwsec.org/scis/2024/program.html#2B4)
<br>
[Experimental analysis of integer factorization methods using lattices](https://www.iwsec.org/2024/program.html)

## NTLライブラリについて
C++のほうでは，NTLライブラリを使用しています．そのため，これを入れていない場合は実行できません．そのため，NTLライブラリを導入していない場合は
```sh
sudo apt-get install -y libntl-dev
```
でインストールしてください．実はインストールされてましたという場合でも
```sh
Reading package lists... Done
Building dependency tree... Done
Reading state information... Done
libntl-dev is already the newest version (11.5.1-1).
0 upgraded, 0 newly installed, 0 to remove and 17 not upgraded.
```
と出るだけで，特にエラーが起こるとかは無いと思うので，導入しているか怪しい方も一応実行しておくとよいかもしれません．

## 実行方法
まず
```sh
make
```
を実行して，a.outを生成します：
```sh
clang++ -Ofast -fopenmp -mtune=native -march=native -pg -std=c++2b lat_fact.cpp -lntl
```

このとき，
```sh
lat_fact.cpp:347:1: warning: non-void function does not return a value in all control paths [-Wreturn-type]
}
^
1 warning generated.
```
というwarningがでることがありますが，こちらについては無視していただいて構いません．

a.outが生成されたら
```sh
./a.out 1 m
```
を実行すると，mビットの合成数の素因数分解が実行されます．以下，その例です
```sh
$ time ./a.out 1 50
1 pairs were found. (800 pairs are needed.)
2 pairs were found. (800 pairs are needed.)
3 pairs were found. (800 pairs are needed.)
4 pairs were found. (800 pairs are needed.)
5 pairs were found. (800 pairs are needed.)
...
535 pairs were found. (800 pairs are needed.)
536 pairs were found. (800 pairs are needed.)
537 pairs were found. (800 pairs are needed.)
X = 36935548444779, Y = 951966701242008
==============================
N = 1122827332192471, bit size = 50
c = 5.0, beta = 2.0
N = 28045841 * 40035431
loop times = 18058
number of sr-pairs = 537
==============================
#Vector = 40648

real    0m22.054s
user    0m21.906s
sys     0m0.130s
```

また，ビット長ではなく，合成数Nを直接指定することもできます：
```sh
./a.out 0 N
```
以下，その実行結果
```sh
time ./a.out 0 1097137578458879
1 pairs were found. (800 pairs are needed.)
2 pairs were found. (800 pairs are needed.)
3 pairs were found. (800 pairs are needed.)
4 pairs were found. (800 pairs are needed.)
5 pairs were found. (800 pairs are needed.)
6 pairs were found. (800 pairs are needed.)
7 pairs were found. (800 pairs are needed.)
...
534 pairs were found. (800 pairs are needed.)
535 pairs were found. (800 pairs are needed.)
536 pairs were found. (800 pairs are needed.)
X = 555427592894032, Y = 788293159767369
==============================
N = 1097137578458879, bit size = 50
c = 5.0, beta = 2.0
N = 55703429 * 19696051
loop times = 17580
number of sr-pairs = 536
==============================
#Vector = 39619

real    0m22.330s
user    0m22.212s
sys     0m0.110s
```
