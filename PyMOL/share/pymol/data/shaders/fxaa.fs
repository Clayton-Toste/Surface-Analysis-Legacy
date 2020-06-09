#version 120

// Based on NVidia FXAA shader

uniform sampler2D color_texture;
uniform vec2 inv_dimensions;

varying vec2 textureLookup ;

#define int2 ivec2
#define float2 vec2
#define float3 vec3
#define float4 vec4
#define FxaaBool3 bvec3
#define FxaaInt2 ivec2
#define FxaaFloat2 vec2
#define FxaaFloat3 vec3
#define FxaaFloat4 vec4
#define FxaaBool2Float(a) mix(0.0, 1.0, (a))
#define FxaaPow3(x, y) pow(x, y)
#define FxaaSel3(f, t, b) mix((f), (t), (b))

#define FXAA_EDGE_THRESHOLD      (1.0/8.0)
#define FXAA_EDGE_THRESHOLD_MIN  (1.0/24.0)
#define FXAA_SEARCH_STEPS        8
#define FXAA_SEARCH_THRESHOLD    (1.0/4.0)
#define FXAA_SUBPIX_CAP          (3.0/4.0)
#define FXAA_SUBPIX_TRIM         (1.0/4.0)

#define FXAA_SUBPIX_TRIM_SCALE (1.0/(1.0 - FXAA_SUBPIX_TRIM))

// Return the luma, the estimation of luminance from rgb inputs.
// This approximates luma using one FMA instruction,
// skipping normalization and tossing out blue.
// FxaaLuma() will range 0.0 to 2.963210702.
float FxaaLuma(float4 rgb) 
{
    return rgb.y * (0.587/0.299) + rgb.x; 
} 

float4 FxaaPixelShader(float2 pos, sampler2D tex, float2 rcpFrame) 
{
/*----------------------------------------------------------------------------
            EARLY EXIT IF LOCAL CONTRAST BELOW EDGE DETECT LIMIT
---------------------------------------------------------------------------*/
    float4 rgbM = texture2D(color_texture, pos.xy );

    float4 rgbN = texture2D(color_texture, pos.xy + rcpFrame * vec2( 0.0,-1.0));
    float4 rgbW = texture2D(color_texture, pos.xy + rcpFrame * vec2(-1.0, 0.0));
    float4 rgbE = texture2D(color_texture, pos.xy + rcpFrame * vec2( 1.0, 0.0));
    float4 rgbS = texture2D(color_texture, pos.xy + rcpFrame * vec2( 0.0, 1.0));

    // If partly transparent, move the nearby pixels towards the middle
    rgbN.rgb = mix( rgbN.rgb, rgbM.rgb, 1. - rgbN.a );
    rgbW.rgb = mix( rgbW.rgb, rgbM.rgb, 1. - rgbW.a );
    rgbE.rgb = mix( rgbE.rgb, rgbM.rgb, 1. - rgbE.a );
    rgbS.rgb = mix( rgbS.rgb, rgbM.rgb, 1. - rgbS.a );

    float lumaN = FxaaLuma(rgbN);
    float lumaW = FxaaLuma(rgbW);
    float lumaM = FxaaLuma(rgbM);
    float lumaE = FxaaLuma(rgbE);
    float lumaS = FxaaLuma(rgbS);

    float rangeMin = min(lumaM, min(min(lumaN, lumaW), min(lumaS, lumaE)));
    float rangeMax = max(lumaM, max(max(lumaN, lumaW), max(lumaS, lumaE)));

    float range = rangeMax - rangeMin;

    float test_range = max(FXAA_EDGE_THRESHOLD_MIN, rangeMax * FXAA_EDGE_THRESHOLD);

    if(range < test_range) {
        return rgbM; 
    }

    float3 rgbL = ( rgbN + rgbW + rgbM + rgbE + rgbS ).rgb;

/*----------------------------------------------------------------------------
                               COMPUTE LOWPASS
---------------------------------------------------------------------------*/

    float lumaL = (lumaN + lumaW + lumaE + lumaS) * 0.25;
    float rangeL = abs(lumaL - lumaM);
    float blendL = max(0.0, 
        (rangeL / range) - FXAA_SUBPIX_TRIM) * FXAA_SUBPIX_TRIM_SCALE; 
    blendL = min(FXAA_SUBPIX_CAP, blendL);
   
/*----------------------------------------------------------------------------
                    CHOOSE VERTICAL OR HORIZONTAL SEARCH
---------------------------------------------------------------------------*/

    float4 rgbNW = texture2D(color_texture, pos.xy + rcpFrame * vec2(-1.0,-1.0));
    float4 rgbNE = texture2D(color_texture, pos.xy + rcpFrame * vec2( 1.0,-1.0));
    float4 rgbSW = texture2D(color_texture, pos.xy + rcpFrame * vec2(-1.0, 1.0));
    float4 rgbSE = texture2D(color_texture, pos.xy + rcpFrame * vec2( 1.0, 1.0));

    // If partly transparent, move nearby pixels towards the middle
    rgbNW.rgb = mix( rgbNW.rgb, rgbM.rgb, 1. - rgbNW.a );
    rgbNE.rgb = mix( rgbNE.rgb, rgbM.rgb, 1. - rgbNE.a );
    rgbSW.rgb = mix( rgbSW.rgb, rgbM.rgb, 1. - rgbSW.a );
    rgbSE.rgb = mix( rgbSE.rgb, rgbM.rgb, 1. - rgbSE.a );

    rgbL += (rgbNW + rgbNE + rgbSW + rgbSE).rgb;
   rgbL *= vec3(1.0/9.0);
    float lumaNW = FxaaLuma(rgbNW);
    float lumaNE = FxaaLuma(rgbNE);
    float lumaSW = FxaaLuma(rgbSW);
    float lumaSE = FxaaLuma(rgbSE);
    float edgeVert = 
        abs((0.25 * lumaNW) + (-0.5 * lumaN) + (0.25 * lumaNE)) +
        abs((0.50 * lumaW ) + (-1.0 * lumaM) + (0.50 * lumaE )) +
        abs((0.25 * lumaSW) + (-0.5 * lumaS) + (0.25 * lumaSE));
    float edgeHorz = 
        abs((0.25 * lumaNW) + (-0.5 * lumaW) + (0.25 * lumaSW)) +
        abs((0.50 * lumaN ) + (-1.0 * lumaM) + (0.50 * lumaS )) +
        abs((0.25 * lumaNE) + (-0.5 * lumaE) + (0.25 * lumaSE));
    bool horzSpan = edgeHorz >= edgeVert;
    float lengthSign = horzSpan ? -rcpFrame.y : -rcpFrame.x;
    if(!horzSpan) lumaN = lumaW;
    if(!horzSpan) lumaS = lumaE;
    float gradientN = abs(lumaN - lumaM);
    float gradientS = abs(lumaS - lumaM);
    lumaN = (lumaN + lumaM) * 0.5;
    lumaS = (lumaS + lumaM) * 0.5;

/*----------------------------------------------------------------------------
                CHOOSE SIDE OF PIXEL WHERE GRADIENT IS HIGHEST
---------------------------------------------------------------------------*/    

    bool pairN = gradientN >= gradientS;
    if(!pairN) lumaN = lumaS;
    if(!pairN) gradientN = gradientS;
    if(!pairN) lengthSign *= -1.0;
    float2 posN;
    posN.x = pos.x + (horzSpan ? 0.0 : lengthSign * 0.5);
    posN.y = pos.y + (horzSpan ? lengthSign * 0.5 : 0.0);
    
/*----------------------------------------------------------------------------
                         CHOOSE SEARCH LIMITING VALUES
----------------------------------------------------------------------------*/

    gradientN *= FXAA_SEARCH_THRESHOLD;
    
/*----------------------------------------------------------------------------
    SEARCH IN BOTH DIRECTIONS UNTIL FIND LUMA PAIR AVERAGE IS OUT OF RANGE
---------------------------------------------------------------------------*/

    float2 posP = posN;
    float2 offNP = horzSpan ? 
        FxaaFloat2(rcpFrame.x, 0.0) :
        FxaaFloat2(0.0, rcpFrame.y); 

    float lumaEndN = lumaN;
    float lumaEndP = lumaN;
    bool doneN = false;
    bool doneP = false;
    posN += offNP * FxaaFloat2(-1.0, -1.0);
    posP += offNP * FxaaFloat2( 1.0,  1.0);
    for(int i = 0; i < FXAA_SEARCH_STEPS; i++) {
        if(!doneN) lumaEndN = FxaaLuma(texture2D(color_texture, posN.xy));
        if(!doneP) lumaEndP = FxaaLuma(texture2D(color_texture, posP.xy));
        doneN = doneN || (abs(lumaEndN - lumaN) >= gradientN);
        doneP = doneP || (abs(lumaEndP - lumaN) >= gradientN);
        if(doneN && doneP) break;
        if(!doneN) posN -= offNP;
        if(!doneP) posP += offNP; 
    }
    
/*----------------------------------------------------------------------------
               HANDLE IF CENTER IS ON POSITIVE OR NEGATIVE SIDE 
---------------------------------------------------------------------------*/

    float dstN = horzSpan ? pos.x - posN.x : pos.y - posN.y;
    float dstP = horzSpan ? posP.x - pos.x : posP.y - pos.y;
    bool directionN = dstN < dstP;
    lumaEndN = directionN ? lumaEndN : lumaEndP;
/*----------------------------------------------------------------------------
         CHECK IF PIXEL IS IN SECTION OF SPAN WHICH GETS NO FILTERING
---------------------------------------------------------------------------*/

    if(((lumaM - lumaN) < 0.0) == ((lumaEndN - lumaN) < 0.0)) 
        lengthSign = 0.0;

/*----------------------------------------------------------------------------
                COMPUTE SUB-PIXEL OFFSET AND FILTER SPAN
---------------------------------------------------------------------------*/

    float spanLength = (dstP + dstN);
    dstN = directionN ? dstN : dstP;
    float subPixelOffset = (0.5 + (dstN * (-1.0/spanLength))) * lengthSign;
    float3 rgbF = texture2D(color_texture, FxaaFloat2(
        pos.x + (horzSpan ? 0.0 : subPixelOffset),
        pos.y + (horzSpan ? subPixelOffset : 0.0))).xyz;
    
    return vec4( mix( rgbF, rgbL, blendL ), rgbM.a );
}

void main()
{
  //  gl_FragColor = vec4(texture2D(color_texture, textureLookup.st).a);
  gl_FragColor = FxaaPixelShader(textureLookup.st, color_texture, inv_dimensions);
	//        FxaaPixelShader(gl_TexCoord[0].st, color_texture, inv_dimensions), 1.0);
}

