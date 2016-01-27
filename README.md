# FLIP_FLUID
This is a flip fluid implement acooding to 3D Particle in Cell / Fluid Implicit Particle Fluid
Solver using OpenMP directives.

Platform:           Windows 10

Develop Tool:       Visual Studio 2013

Program language:   C++  OpenGL

External libraries: glfw glew

usage:

1.Download zip or clone this repository.

2.Open "FlipFluid.sln" using Visual Studio 2013(Only support version 2013 now).

3.Change the Configuration to Release.

4.Run.
In this step you may see "glew32.dll are not in your compute" something like this,
go "lib" sub folder and copy the "glew32.dll" file to "Release" folder.After this,
it should work.
  
You can umcomment "#define (someeffect)" to see different effect.
For surface reconstrution I use the code form https://code.google.com/archive/p/flip3d/.