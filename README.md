
# ell_ell_isogeny_sage

今回の論文の疑似コードに沿ったコードです. 現在, 結果が整合的でないです. 

## 使い方

最初にload_this.sage内のimportの部分とloadの部分をコピペしてください. 

1. B-SIDH attackを実装する場合は, 同ファイルの対応する部分(#attack for B-SIDH)をコピペして下さい. (途中でassertion errorを起こします: 点像のorderが正しくないとなります.)

2. 例を計算する場合は, 同ファイルの対応する部分(#calculation example)をコピペしてください. (途中でassertion errorを起こします: 点像のorderが正しくないとなります.)
   例は定義域が楕円曲線の直積で, 素数次数の同種写像を計算しています.

3. countはまだ整備されていません. 

## コードの説明

class_count.pyがcountするための体を作るコードが書かれています. 

calss_theta.pyはtheta座標について書かれています. テータ座標の演算などもここです. 

func_elliptic.pyは楕円曲線に関するアルゴリズムです. ほとんどもう一つのBSIDH_attack_sageと同じです. 

func_for_attack.pyはB-SIDH attackに関することが書かれています. 

func_fraction.pyは分数の計算について書かれています. 分数をよく扱うので, それ関連のものです. 

func_isogeny.pyは同種写像の計算について書かれています. これがメインです. 

# BSIDH_attack_sage   

同じくsageのコードで正しく動きますが, 割り算を用いています. 

## 使い方

Read_me.pyの一番上の`load("setting.py")`をコピペしてください. 

1. B-SIDH攻撃したい場合は下の部分をコピペしてください.

2. 計算時間を測りたい場合はその下の部分をコピペしくてださい. 




