#version 120

attribute vec3 a_Vertex;
varying vec2 textureLookup ;

void main(void)
{
  gl_Position = vec4(a_Vertex.x, a_Vertex.y, .5, 1.);
  textureLookup = (1. + a_Vertex.xy) / 2.;
}

