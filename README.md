# lattice_factorization

このリポジトリに含まれているソースたちは，大学院での研究の際に作成されたものです．以下，その論文：
<br>[
近似最近ベクトル探索と埋め込み法を用いた格子による素因数分解法の実装報告](https://www.iwsec.org/scis/2024/program.html#2B4)

[Experimental analysis of integer factorization methods using lattices](https://link.springer.com/chapter/10.1007/978-981-97-7737-2_8)

## 使用方法

例として，Pythonで使用する``lat_fact.py``及び``test.py``が同梱されているので，そちらを実行した際の挙動を見てみます．本ファイルを実行すると,以下のような挙動になるかと思います．（Windows，Linuxどちらでも同様です）

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

test.pyでは，格子を用いた素因数分解法を用いて，ランダムに生成された50ビットRSA型合成数の素因数分解を行っています．最終的に一番下に出てきた整数の組が素因数の組になっています．その上に表示されている謎の文字列たちは格子素因数分解をするにあたってどの程度進捗しているかなどを表している，いわば格子素因数分解法の情報です．``lat_fact.py``内で定義されている``lat_fact``は，第一引数に素因数分解したい合成数，もしくはそのビットサイズを，第二引数にビットサイズで入力するか，合成数をそのまま入力するか，第三引数に進捗情報を出力するかをそれぞれ入力します．初期値では，ビットサイズで入力し，進捗情報は出力するようになっています．

続いて，dllファイルもしくはsoファイルについて説明します．これは，同梱されている``lat_fact.hpp``及び``lat_fact.cpp``をコンパイルしたものになります．ここでは，``lattice_factorization``関数のみが使えるようになっており，第一引数にビットで入力するか否かを，第二引数に素因数分解したい合成数，もしくはそのビットサイズ（第一引数によって変わります）を，第三引数に情報を出力するか否かを入力します．他の関数を利用したい場合は，適宜C++コードを変更して使用可能にしていただいても構いません．


# 動作確認環境
- Windows 11
- Ubuntu 22.04.5 LTS
