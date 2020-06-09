uniform sampler3D volumeTex;
uniform sampler2D colorTex2D;
uniform sampler1D colorTex1D;
uniform float volumeScale;
uniform float volumeBias;
uniform float sliceDist;
uniform mat4 TexMatrix;
uniform vec3 eyeposM;
uniform vec3 vDirM;
uniform sampler3D carvemask;
uniform bool carvemaskFlag;

varying vec3 vertexM;
varying float fog;
varying vec2 bgTextureLookup;

#include anaglyph_header.fs
#include compute_fog_color.fs

bool iscarvemasked(vec3 t) {
  return carvemaskFlag && texture3D(carvemask, t).r > 0.5;
}

void main()
{
#ifdef volume_mode
#ifdef ortho
  vec3 tex0 = gl_TexCoord[0].xyz;
  vec3 tex1 = gl_TexCoord[1].xyz;
#else // ortho

  // normalized eye to vertex vector in model space
  vec3 eyeToVert = normalize(vertexM - eyeposM);

  // vector from vertex to the back slab plane
  vec3 slabOffset = eyeToVert * (sliceDist / dot(vDirM, eyeToVert));

  // back slab and front slab points in model space
  vec3 sB = vertexM - slabOffset;
  vec3 sF = vertexM + slabOffset;

  // texture coordinates
  vec3 tex0 = vec3(vec4(sB, 1.) * TexMatrix);
  vec3 tex1 = vec3(vec4(sF, 1.) * TexMatrix);

#endif // ortho

  if (iscarvemasked(tex0))
    discard;

  // texture lookup of map values
  vec2 v = vec2(
    texture3D(volumeTex, tex0).r,
    texture3D(volumeTex, tex1).r
  );

  // color lookup in (preintegrated) table
  v = v * volumeScale + volumeBias;
  vec4 color = texture2D(colorTex2D, v);

#else // volume_mode

  if (iscarvemasked(gl_TexCoord[0].xyz))
    discard;

  float v = texture3D(volumeTex, gl_TexCoord[0].xyz).r;
  v = v * volumeScale + volumeBias;
  if (v < 0. || v > 1.) discard;
  vec4 color = texture1D(colorTex1D, v);
#endif // volume_mode

  if (color.a == 0.0)
    discard;

  color = ApplyColorEffects(color, gl_FragCoord.z);

  gl_FragColor = ApplyFog(color, fog);

  PostLightingEffects();
}

