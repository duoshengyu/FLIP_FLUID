# If the text is garbled, please convert encoding to UTF-8.
# FLIP 流体
このプログラムは[3D Particle in Cell / Fluid Implicit Particle Fluid
Solver using OpenMP directives]によって作成したものです。

「SIGGRAPH 2007 course notes FLUID SIMULATION」 共役勾配法の実装


実行ファイル：Release/FlipFluid.exe

*実行環境等を含めた実行方法:

動作確認：Windows 10　Visual Studio 2013

開発時間と人数:　四か月前  二週　一人

開発言語：C++  OpenGL

外部ライブラリ：glfw glew

実行方法："FlipFluid.sln"　Visual Studio 2013で開けたら、F5を押す。

操作方法：

左マウスクリック+ドラッグ　　　視点を回転

ボタン'1'を押すと　　　　　　シン　one dambreak 再現
ボタン'2'を押すと　　　　　　シン　two dambreak 再現


*プログラムを作成する上で苦労した箇所は？
デバッグ、FLIP流体体シミュレーションの勉強

*力をいれて作った部分で、「プログラム上」で特に注意してみてもらいたい箇所は？
粒子とグリッド間のスピードマッピング、FLIP.cppのP2G、G2P 関数
プレッシャーの連立一次方程式を解くための共役勾配法の実装
FLIP.cppのsolvePressureAndNewV()　と　solver.h solver.cpp クラス(flip3dのフローワークによって修正)


*参考にしたソースファイルがあるなら、どの様なところを参考にしましたか？またその部分のファイル名を書いてください
flip3d：https://code.google.com/archive/p/flip3d/
メッシュの構成(surface construction)に参考。

*他人のコードと関数：FLIP_other.cpp　粒子の後処理

                  CIsoSurface.h CIsoSurface.cpp implict.h 
                  Vectors.h Vectors.cpp Sorter.h Sorter.cpp
                  メッシュの構成(surface construction)に使用、
                  提供している実行ファイル、またはそのままプロジェクトをビルドすると
                  使っていません。
───
このプロジェクトはhttps://github.com/duoshengyu/FLIP_FLUIDにも見えます。