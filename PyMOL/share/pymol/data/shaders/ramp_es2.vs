attribute vec4 a_Vertex;
attribute vec4 a_Color;
attribute vec3 a_Normal;

varying vec4 COLOR ;
varying vec3 NORMAL ;

uniform vec3 offsetPt;

uniform mat4 g_ProjectionMatrix;
uniform mat4 g_ModelViewMatrix;

void main()
{
  COLOR = a_Color;
  NORMAL = a_Normal;
  vec4 vertex = a_Vertex + vec4(offsetPt,0.);
  gl_Position = g_ProjectionMatrix * g_ModelViewMatrix * vertex;
}
