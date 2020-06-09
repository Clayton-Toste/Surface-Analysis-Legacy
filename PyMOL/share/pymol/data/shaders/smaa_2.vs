#version 120
#extension GL_EXT_gpu_shader4 : require

uniform vec4 SMAA_RT_METRICS;
uniform float isOutput;

attribute vec3 a_Vertex;

varying vec2 texcoordAttr;
varying vec2 pixcoord;
varying vec4 offset[3];

#include smaa_gen.vs

void main()
{
	vec2 vert = vec2((1. + a_Vertex.xy) / 2.);
	texcoordAttr = vert.xy;
	SMAABlendingWeightCalculationVS(texcoordAttr, pixcoord, offset);
	gl_Position = vec4(a_Vertex.x, a_Vertex.y, 0., 1.);
}
