このコードは、Twitter
https://twitter.com/madnoda/status/1366886071288492034
の一連の流れで、倍精度浮動小数点計算と四倍精度浮動小数点計算の計算結果の違いを示すためのものである

プログラムは、Ubuntu 18.04上で動作確認したが、おそらく、gcc(g++)とgnuplotがインスールされていれば、どんなLINUXや他のOSでも動作すると思われる

$ make -f solar_earth_L2_002.mak 
で、倍精度浮動小数点計算

$ make -f solar_earth_L2_002q.mak 
で、四倍精度浮動小数点計算
を行う

両方とも、gnuplotによるグラフが最終的に作成される

参考までに Ryzen 5 3600 / Ubuntu 18.04 での実行時間は、倍精度が17分、四倍精度が24分であった

なお、これらのコードは、MITライセンスとするが、まだ未完成であり、現時点では、あまり拡散しないでもらいたい

また、コードの不具合等での損害等に関しては、免責とさせてもらう

2021.03.04 野田篤司 Twitter: @madnoda

