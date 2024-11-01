# lattice_factorization

このリポジトリに含まれているソースたちは，大学院での研究の際に作成されたものです．以下，その論文：
<br>[
近似最近ベクトル探索と埋め込み法を用いた格子による素因数分解法の実装報告](https://www.iwsec.org/scis/2024/program.html#2B4)
<br>
[Experimental analysis of integer factorization methods using lattices](https://link.springer.com/chapter/10.1007/978-981-97-7737-2_8)

## 使用方法

基本的にはPythonで使用することを想定しているので，その方法で説明していきます．
まず，例として同梱されている``test.py``が何をしているのかについて説明します．
本ファイルを実行すると,以下のような挙動になるかと思います．

```sh
$ python3 test.py
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
(28045841, 40035431)
```

test.pyでは，格子を用いた素因数分解法を用いて，ランダムに生成された50ビットRSA型合成数の素因数分解を行っています．
最終的に一番下に出てきた整数の組が素因数の組になっています．
その上に表示されている謎の文字列たちは格子素因数分解をするにあたってどの程度進捗しているかなどを表している，いわば格子素因数分解法の情報です．

testを実行するのみではなく，自分で色々やったりしてみたいという方もいるかと思います．そのために``lat_fact.py``内の``lat_fact``関数について説明をしておきます．
当該関数では引数が``N``，``bit_input``，``print_info``の3つがあります．
``N``は素因数分解したい整数，もしくはそのビットサイズです．``bit_input``として1を指定した場合はビットサイズでの入力，0を指定した場合は整数での入力となります．
``print_info``は``* pairs were found``とかそのへんの情報を出力するかを指定するもので，1を指定した場合は情報を出力，0を指定した場合は情報を出力しない様になります．

## 正常に動作しない場合
まず，前提としてWindows上での動作は想定していませんので，Windows上では正常に動作しないと思って下さい．ubuntuやlinux上での動作を想定しています．
その上で，正常に動かない場合は以下を試してみて下さい．

### NTLライブライのインストール

C++のほうでは，NTLライブラリを使用しています．そのため，これを入れてない場合動かないことがあるかもしれません．そのため，NTLライブラリを導入していなくて正常に動作しない場合は
```sh
$ sudo apt-get install -y libntl-dev
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

### soファイルを作る
soファイルをご自身でコンパイルして作っていただくことによって解決するかもしれません．この場合，自分で``g++ hoge``などとコマンドを打たずとも，Makefileを実行いただくことでコンパイルできるようになっているので，
```sh
$ make
```
を実行して，libfact.soを生成して下さい．
```sh
g++ -shared -fPIC -O3 -mavx2 -fopenmp -mtune=native -march=native -mfpmath=both -unroll-loops -o libfact.so src/lat_fact.cpp -lntl
```

また，コンパイルには必ずNTLライブラリが必要になりますので，インストールは済ませておいて下さい．

このとき，
```sh
src/lat_fact.cpp: In function ‘std::vector<long long int> ENUM_CVP(std::vector<std::vector<double> >, std::vector<double>, double, std::vector<double>)’:
src/lat_fact.cpp:409:1: warning: control reaches end of non-void function [-Wreturn-type]
  409 | }
      | ^
```

というwarningがでることがありますが，こちらについては無視していただいて構いません．