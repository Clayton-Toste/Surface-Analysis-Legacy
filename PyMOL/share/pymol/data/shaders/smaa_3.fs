#version 120
#extension GL_EXT_gpu_shader4 : require

uniform vec4 SMAA_RT_METRICS;
uniform sampler2D colorTex;
uniform sampler2D blendTex;
uniform float isOutput;

varying vec2 texcoordAttr;
varying vec4 offset;

#include smaa_gen.fs

void main()
{
	gl_FragColor = SMAANeighborhoodBlendingPS(texcoordAttr, offset, colorTex, blendTex);
}
