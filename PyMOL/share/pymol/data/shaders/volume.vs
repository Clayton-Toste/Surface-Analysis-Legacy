#include webgl_header.vs

uniform vec3 slabTexOffset;

varying vec3 vertexM;

varying float fog;
varying vec2 bgTextureLookup;

void main()
{
  vec4 vertex = g_ModelViewMatrix * gl_Vertex;
  gl_TexCoord[0] = gl_MultiTexCoord0;

#ifdef volume_mode
#ifdef ortho
  gl_TexCoord[1] = gl_TexCoord[0];
  gl_TexCoord[0].zyx -= slabTexOffset;
  gl_TexCoord[1].zyx += slabTexOffset;
#else // ortho

  // vertex varying in model space
  vertexM = vec3(gl_Vertex);

#endif // ortho
#endif // volume_mode

  gl_Position = ftransform();
  fog = (g_Fog_end + vertex.z) * g_Fog_scale;
  bgTextureLookup = (gl_Position.xy/gl_Position.w) / 2.0 + 0.5;
}
