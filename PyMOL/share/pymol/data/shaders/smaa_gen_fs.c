#define texture texture2D
#define textureLod texture2DLod
#define textureLodOffset texture2DLodOffset

#define out inout
#define SMAA_INCLUDE_PS 1
#define SMAA_GLSL_3 1

uniform float SMAA_THRESHOLD;
uniform int SMAA_MAX_SEARCH_STEPS;
uniform int SMAA_MAX_SEARCH_STEPS_DIAG;
uniform int SMAA_CORNER_ROUNDING;

//#define SMAA_PRESET_LOW 1
//#define SMAA_PRESET_HIGH 1
//#define SMAA_PRESET_ULTRA 1

/*
// This is what MAESTRO uses
#define SMAA_INCLUDE_PS 1
#define SMAA_CUSTOM_SL
#define SMAA_PRESET_LOW
#define SMAATexture2D(tex) sampler2D tex
#define SMAATexturePass2D(tex) tex
#define SMAASampleLevelZero(tex, coord) texture2DLod(tex, coord, 0.0)
#define SMAASampleLevelZeroPoint(tex, coord) texture2DLod(tex, coord, 0.0)
#define SMAASampleLevelZeroOffset(tex, coord, offset) texture2DLodOffset(tex, coord, 0.0, offset)
#define SMAASample(tex, coord) texture2D(tex, coord)
#define SMAASamplePoint(tex, coord) texture2D(tex, coord)
#define SMAASampleOffset(tex, coord, offset) texture2D(tex, coord, offset)
#define SMAA_FLATTEN
#define SMAA_BRANCH

#define lerp(a, b, t) mix(a, b, t)
#define saturate(a) clamp(a, 0.0, 1.0)
#define mad(a, b, c) (a * b + c)

#define float2 vec2
#define float3 vec3
#define float4 vec4
#define int2 ivec2
#define int3 ivec3
#define int4 ivec4
#define bool2 bvec2
#define bool3 bvec3
#define bool4 bvec4

*/
#include "SMAA.hlsl"

