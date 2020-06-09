#version 120
#extension GL_EXT_gpu_shader4 : require

uniform vec4 SMAA_RT_METRICS;

uniform sampler2D edgesTex;
uniform sampler2D areaTex;
uniform sampler2D searchTex;
uniform float isOutput;

varying vec2 texcoordAttr;
varying vec2 pixcoord;
varying vec4 offset[3];

#include smaa_gen.fs

void main()
{
	gl_FragColor = SMAABlendingWeightCalculationPS(texcoordAttr, pixcoord, offset, edgesTex, areaTex, searchTex, vec4(0,0,0,0));
}
